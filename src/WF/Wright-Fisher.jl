# Definitionn of the Species.
mutable struct Species
    n::Int
    f::Float32
    sequence::Array{Int}
end

make_clone(s::Species) = Species(copy(s.n), copy(s.f), copy(s.sequence));

function get_mutants(s0::Species, mu, L, h)
    newSpecies = []
    
    if s0.n > 0
        # Select the number of mutants
        nMut = rand(Binomial(s0.n, mu * L))
        if nMut > 0
            # Remainig original should be the A1.n - nMut
            s0.n -= nMut 
            # Sites where that involve mutations
            sites = rand(1:L, nMut)

            for i in sites
                s = make_clone(s0)
                s.n = 1
                s.sequence[i] = 1 - s.sequence[i] # flip 0 <-> 1
                s.f           = 1. + sum(h .* s.sequence)
                newSpecies = push!(newSpecies, s)
            end
            # Adding the s0 it selfs if there are remaining s0, that is s0.n>0 .
            if s0.n > 0
                newSpecies = push!(newSpecies, s0)
            end
        end
        if nMut == 0 
            newSpecies = [make_clone(s0)]
        end
    end
    return newSpecies
end

function get_allele_freq(pop, L, N)
    f1 = zeros(L)
    for i in 1:L
        for s in pop    
            f1[i] += s.n * s.sequence[i] / N 
        end
    end
    return f1
end;

function get_allele_freq_pairwise(pop, L, N)
    f1 = zeros(L)
    f2 = zeros(L,L)
    
    for s in pop    
        f1 += s.n * s.sequence[:] / N
        f2 += s.n * s.sequence[:] * s.sequence[:]' / N
    end
    
    return (f1,f2)
end;

function get_allels_freq_series(
        h,
        ensemble_id,
        dir_out, 
        N, 
        sampling_interval_set; 
        file_key="standard",
        mu=1e-3, 
        T=301, 
        n_sampling_max=100)
    
    L = length(h)
    #--------- start WF model ------------#
    tStart = 0
    tEnd = T
    
    pop, sVec, nVec = [], [], [];
    f_ini = 0.001 # it should be >0 
    #pop = [Species(N, f_ini, zeros(L))] #SIMPLE
    pop = [Species(Int(N/5), f_ini, rand(0:1, L))] # COMPLEX
    for i in 1:4
        pop = push!(pop, Species(Int(N/5), f_ini, rand(0:1, L))) # COMPLEX
    end

    sVec = [Array([zeros(L)])]
    nVec = [Array([N])];
    allels_series = [];
    f1 = get_allele_freq(pop, L, N)
    allels_series = copy(f1);

    fname_out = @sprintf("%s./allele-traject_id-%d_mu-%.3f_N-%d_WF_%s.txt",
                            dir_out, ensemble_id, mu, N, file_key)
    fout = open(fname_out, "w")
    
    t_start = time();
    for t in tStart:tEnd
        r = Array([s.n * s.f for s in pop]) * 1.0 
        prob_r = r / sum(r)
        if(length(prob_r)>1)
            n = rand(Multinomial(N, prob_r))
        end  
        if(length(prob_r)==1)
            n = [N]
        end
        newPop = []

        for i in 1:length(pop)
            pop[i].n = n[i]
            p = get_mutants(pop[i], mu, L, h)

            for j in 1:length(p)
                unique = true
                for k in 1:length(newPop)
                    if( p[j].sequence == newPop[k].sequence)
                        unique = false
                        global newPop[k].n += p[j].n
                        break
                    end
	    end
                if unique
                    newPop = push!(newPop, p[j])
                end
            end
        end
        pop = newPop

        f1 = get_allele_freq(pop, L, N)
        allels_series = hcat(allels_series, f1)
        # Output the trajectory to a file.  
        if(t in sampling_interval_set) 
           for s in pop
                lineout = @sprintf "%d\t%d\t%.4f\t%s" t s.n s.f join([x for x in s.sequence], " ")
                println(fout, lineout)
            end
        end
    end
    t_end = time()
    close(fout)
    return allels_series
end;

"""
T = 301 # stopping time. 
K = 6 # number of observation including begning and end. 
"""
function get_random_sampling_interval(T=300, K=6)    
    interval_raw = rand(Uniform(0,1), K)
    normalize_interval = sum(interval_raw)
    dt_list = [Int(floor(T * interval_raw[i]/normalize_interval)) for i in 1:(K-1)]
    dt_last = T-sum(dt_list)
    push!(dt_list, dt_last);
    
    dt_accum_list = []
    dt_accum = 1
    push!(dt_accum_list, dt_accum)
    
    for i in 1:K
        dt_accum += dt_list[i]
        push!(dt_accum_list, dt_accum)
    end
    return (copy(dt_list), copy(dt_accum_list))
end;
