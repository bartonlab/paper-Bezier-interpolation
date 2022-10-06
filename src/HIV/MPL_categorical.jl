using Pkg
using DelimitedFiles
using Profile    
using Random
using LinearAlgebra
using Plots
using Printf
rng = Random.MersenneTwister(1234);

km(i,a,q) = (i-1)*q+a

function get_x1_x2(q, L, n_t, sample_t)
    Neff = Int(sum(n_t));
    n_species = size(sample_t,2);
    x1 = zeros(L*q);
    x2 = zeros(L*q,L*q);
    for m in 1:n_species
        x1 += n_t[m] / Neff * sample_t[:,m]
        x2 += n_t[m] / Neff * sample_t[:,m] * sample_t[:,m]'
    end 
    return (x1, x2)
end
 
function get_integrated_corr_mean(t_local, q, L, x1_t1, x1_t2, x2_t1, x2_t2)
    x1_out = 0.5 * (x1_t1 + x1_t2)
    
    C = 0.5 * (x2_t1 + x2_t2) 
    C += - (1.0/3) * (x1_t1 * x1_t1' + x1_t2 * x1_t2' )
    mat_temp = x1_t1 * x1_t2'
    C += - (1.0/6) * (mat_temp + mat_temp' )
    
    C_diag = 0.5 * (x1_t1 + x1_t2) - (1.0/3) * ( x1_t1.^2 + x1_t1 .* x1_t2 + x1_t2.^2 )
    C[diagind(C)] = C_diag
    
    return (x1_out, C)
end;

function get_flux(q, L, x1, muMat)
	fluxIn = zeros(q*L)
	fluxOut = zeros(q*L)
	for i in 1:L
		fluxIn[km.(i,1:q,q)] = muMat' * x1[km.(i,1:q,q)]
		fluxOut[km.(i,1:q,q)] = x1[km.(i,1:q,q)] .* ( muMat * ones(q) ) 
	end
	return (fluxIn, fluxOut)
end;


function get_flux_poly(qL, q, L, x1, muMat, 
		       polymorphic_positions, 
		       polymorphic_states_set,
		       original_index_to_polymorophic_index) 
	fluxIn = zeros(qL)
	fluxOut = zeros(qL)
	len_poly_positions = size(polymorphic_positions,1)
	for n in 1:len_poly_positions  
		i = polymorphic_positions[n]
		for a in polymorphic_states_set[n]
			for b in polymorphic_states_set[n]
				# effective index using mapping: (i,a) -> u
				u = original_index_to_polymorophic_index[km(i,a,q)]
				v = original_index_to_polymorophic_index[km(i,b,q)]
				@assert(u>0 && v>0)
				fluxIn[u] += muMat[b,a] * x1[v]
				fluxOut[u] += x1[u] *  muMat[a,b] 
			end	
		end
	end
	return (fluxIn, fluxOut)
end;


function get_x1_x2_psa(qL, n_t, sample_t, scaling_psc, pseud_count=0.01)
	alpha = 1-pseud_count
	Neff = Int(sum(n_t));
	n_species = size(sample_t,2);
	x1 = zeros(qL);
	x2 = zeros(qL,qL);
	for m in 1:n_species
		x1 += n_t[m] / Neff * sample_t[:,m]
		x2 += n_t[m] / Neff * (sample_t[:,m]*sample_t[:,m]')
	end
	ones_m = ones(qL)
	pseudo_count_vector = ((1-alpha)*scaling_psc) .* ones_m
	x1 = alpha * copy(x1) + pseudo_count_vector
	pseudo_count_mat = pseudo_count_vector * pseudo_count_vector'
	pseudo_count_mat[diagind(pseudo_count_mat)] = pseudo_count_vector
	x2 = alpha * copy(x2) + pseudo_count_mat    
	return (x1, x2)
end


"""
function get_x1_x2_psa(q, L, n_t, sample_t, pseud_count=2)
    Neff = Int(sum(n_t));
    pseud_count /= Neff 
    pseud_count = minimum([pseud_count, 0.3]) 
    @printf( "effective pseud_count = %.4f\n", pseud_count) 
    alpha = 1-pseud_count
    n_species = size(sample_t,2);
    x1 = zeros(L*q);
    x2 = zeros(L*q,L*q);
    for m in 1:n_species
        ones_m = ones(size(sample_t[:,m]))
        # the ones term should be scaled by the number of possible states.
	data_in = copy(alpha*sample_t[:,m] + (1-alpha)*ones_m)
        x1 += n_t[m] / Neff * data_in
        x2 += n_t[m] / Neff * (alpha*sample_t[:,m]*sample_t[:,m]' + (1-alpha)*ones_m*ones_m')
    end 
    return (x1, x2)
end;
"""

# Get integrated covariance matrix and drift. 
function get_integrated_Ctot_drift_x1mean_simple(q, L, data, time_list_in, muMat, pseud_count)
    time_list = unique(time_list_in) 
    C_tot = zeros(L*q,L*q)
    drift = zeros(L*q)
    x1_mean = zeros(L*q)
    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1_t1, x2_t1) = get_x1_x2_psa(q*L, n_t, sample_t, pseud_count)
    #(x1_t1, x2_t1) = get_x1_x2(q, L, n_t, sample_t)
    x1_init = copy(x1_t1)
    x1_mean = copy(x1_t1) * t_old
    Dt_sum = 0 

    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t); # without interpolate
        (x1_t2, x2_t2) = get_x1_x2_psa(q*L, n_t, sample_t, pseud_count)
        #(x1_t2, x2_t2) = get_x1_x2(q, L, n_t, sample_t)
        (x1_integrated, C_integreated) = get_integrated_corr_mean(t, q, L, x1_t1, x1_t2, x2_t1, x2_t2)
        x1_t1 = copy(x1_t2); 
        x2_t1 = copy(x2_t2)

        Dt = t - t_old; 
        t_old = t
        (fluxIn, fluxOut) = get_flux(q, L, x1_integrated, muMat) 
        drift += Dt * (fluxIn - fluxOut)
        C_tot += Dt * C_integreated
        x1_mean += Dt * copy(x1_t2) 
        Dt_sum += Dt
    end
    x1_mean /= (Dt_sum+time_list[1]) 
    x1_fini = copy(x1_t1)
    return (C_tot, drift, x1_mean, x1_init, x1_fini)
end;

# Get integrated covariance matrix and drift. 
function get_integrated_Ctot_drift_x1mean_simple_poly(
		poly_seq_len, q, L, data, time_list_in, muMat, pseud_count,
		polymorphic_positions, 
		polymorphic_states_set,
		original_index_to_polymorophic_index,
		scaling_psc
		)
    
    time_list = unique(time_list_in) 
    C_tot = zeros(poly_seq_len, poly_seq_len)
    drift = zeros(poly_seq_len)
    x1_mean = zeros(poly_seq_len)
    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1_t1, x2_t1) = get_x1_x2_psa(poly_seq_len, n_t, sample_t, scaling_psc, pseud_count)
    #(x1_t1, x2_t1) = get_x1_x2(q, L, n_t, sample_t)
    x1_init = copy(x1_t1)
    x1_mean = copy(x1_t1) * t_old
    Dt_sum = 0 
    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t); # without interpolate
        (x1_t2, x2_t2) = get_x1_x2_psa(poly_seq_len, n_t, sample_t, scaling_psc, pseud_count)
        #(x1_t2, x2_t2) = get_x1_x2(q, L, n_t, sample_t)
        (x1_integrated, C_integreated) = get_integrated_corr_mean(t, q, L, x1_t1, x1_t2, x2_t1, x2_t2)
        x1_t1 = copy(x1_t2); 
        x2_t1 = copy(x2_t2)

        Dt = t - t_old; 
        t_old = t
        (fluxIn, fluxOut) = get_flux_poly(poly_seq_len, q, L, x1_integrated, muMat, 
				polymorphic_positions, 
			       	polymorphic_states_set,
			       	original_index_to_polymorophic_index)
        
	drift += Dt * (fluxIn - fluxOut)
        C_tot += Dt * C_integreated
        x1_mean += Dt * copy(x1_t2) 
        Dt_sum += Dt
    end
    x1_mean /= (Dt_sum+time_list[1]) 
    x1_fini = copy(x1_t1)
    return (C_tot, drift, x1_mean, x1_init, x1_fini)
end;


#function selction_coeff_cholfact(mu, drift, C_tot, x1_init, x1_fini)
#    return s = inv(cholesky(C_tot'*C_tot)) * C_tot' * (x1_fini-x1_init-mu*drift)
#end;

function selction_coeff_cholfact(mu, reg, drift, C_tot, x1_init, x1_fini)
    return s = (C_tot+reg*I) \ (x1_fini-x1_init-mu*drift)
end;

function get_MPLseq(fname_in)
    data_raw = readdlm(fname_in)
    specific_times = Int.(data_raw[:,1])
    unique_time = unique(specific_times)
    index_sorted_time = sortperm(specific_times)
    MPLseq_raw = []
    L = size(data_raw,2)-2
    q = size(unique(data_raw[:, 3:end]), 1)
    #@show q # For DNA it should be 5
    N = size(data_raw,1)
    @printf("L=%d, q=%d, N=%d\n", L, q, N)
   
    for n in 1:N
        temp1 = []
        temp1 = vcat(temp1, specific_times[n])
        temp1 = vcat(temp1, data_raw[n,2])
        temp1 = vcat(temp1, 0.9); 
        for i in 1:L
            temp2 = zeros(Int, 5); 
            a = Int(data_raw[n,2+i]) + 1
            temp2[a] = 1
            temp1 = vcat(temp1, copy(temp2))
        end
        if(length(MPLseq_raw)>0)
            MPLseq_raw = hcat(MPLseq_raw, copy(temp1))
        end
        if(length(MPLseq_raw)==0)
            MPLseq_raw = copy(temp1)
        end
    end
    MPLseq_raw = copy(MPLseq_raw');
    return (q,L,N,MPLseq_raw,specific_times)
end;

