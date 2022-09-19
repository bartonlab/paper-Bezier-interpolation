#Pkg.add("Pkg"); import Pkg; 
using Pkg
using Distributed
using DelimitedFiles
using Distances
using StatsBase 
using Profile    
using Random
using Statistics
using LinearAlgebra
using Distributions
using Plots
#using StatsPlots
using Printf
rng = Random.MersenneTwister(1234);

function get_x1_x2(L, n_t, sample_t)
    Neff = Int(sum(n_t))
    x1 = zeros(L)
    x2 = zeros(L,L);
    n_species = size(sample_t,2);
    for m in 1:n_species
        x1 += n_t[m] / Neff * sample_t[:,m]
        x2 += n_t[m] / Neff * sample_t[:,m] * sample_t[:,m]'
    end

    return (x1, x2)
end;

function get_x1_x2_C(L, n_t, sample_t)
    Neff = Int(sum(n_t));
    n_species = size(sample_t,2);
    x1 = zeros(L)
    x2 = zeros(L,L);
    for m in 1:n_species
        x1 += n_t[m] / Neff * sample_t[:,m]
        x2 += n_t[m] / Neff * sample_t[:,m] * sample_t[:,m]'
    end

    C = copy(x2)  - x1 * x1';
    for i in 1:L
        C[i,i] = x1[i]*(1-x1[i])
    end
    return (x1, x2, C)
end;


function get_mean_corr(t_local, L, x1, x2)
    C =  copy(x2 - x1 * x1');
    for i in 1:L
        C[i,i] = x1[i]*(1-x1[i])
    end
    return (x1,C)
end

function get_interpolate_mean_corr(t_local, L, x1_t1, x1_t2, x2_t1, x2_t2)
    r = 0.5
    x1 = r * (x1_t1 + x1_t2) 
    x2 = r * (x2_t1 + x2_t2) 
    C = copy(x2) - x1 * x1';
    
    C_diag = x1 .* (ones(size(x1)) - x1)
    C[diagind(C)] = C_diag
    return (x1, C)
end;

function get_integrated_corr_mean(t_local, L, x1_t1, x1_t2, x2_t1, x2_t2)
    x1_out = 0.5 * (x1_t1 + x1_t2)
    C = 0.5 * (x2_t1 + x2_t2) 
    C += - (1.0/3) * (x1_t1 * x1_t1' + x1_t2 * x1_t2' )
    mat_temp = x1_t1 * x1_t2'
    C += - (1.0/6) * (mat_temp + mat_temp' )
    
    C_diag = 0.5 * (x1_t1 + x1_t2) - (1.0/3) * ( x1_t1.^2 + x1_t1 .* x1_t2 + x1_t2.^2 )
    C[diagind(C)] = C_diag
    
    return (x1_out, C)
end;

function get_x1_x2_C_interpolate(L, n_t1, n_t2, sample_t1, sample_t2)
    Neff_t1 = Int(sum(n_t1))
    Neff_t2 = Int(sum(n_t1))
    x1_t1 = zeros(L); x1_t2 = zeros(L)
    x2_t1 = zeros(L,L); x2_t2 = zeros(L,L);
    
    for m in 1:size(sample_t1,2)
        x1_t1 += n_t1[m] / Neff_t1  * sample_t1[:,m] 
        x2_t1 += n_t1[m] / Neff_t1 * sample_t1[:,m] * sample_t1[:,m]'
    end
    for m in 1:size(sample_t2,2)
        x1_t2 += n_t2[m] / Neff_t2  * sample_t2[:,m] 
        x2_t2 += n_t2[m] / Neff_t2 * sample_t2[:,m] * sample_t2[:,m]'
    end
    r = 0.5 
    x1 = r * (x1_t1 + x1_t2) 
    x2 = r * (x2_t1 + x2_t2) 
    C = copy(x2) - x1 * x1';
    for i in 1:L
        C[i,i] = x1[i]*(1-x1[i])
    end
    return (x1, x2, C)
end;
#
#
## Integration should be done for both the flux and covariance terms. 
# This function does not yet treat the integration of the flux. 
function get_x1_x2_C_integrate(L, n_t1, n_t2, sample_t1, sample_t2)
    Neff_t1 = Int(sum(n_t1));
    Neff_t2 = Int(sum(n_t2));
    x1_t1 = zeros(L); x1_t2 = zeros(L)
    x2_t1 = zeros(L,L); x2_t2 = zeros(L,L);
    for m in 1:size(sample_t1,2)
        x1_t1 += n_t1[m] / Neff_t1  * sample_t1[:,m] 
        x2_t1 += n_t1[m] / Neff_t1 * sample_t1[:,m] * sample_t1[:,m]'
    end
    for m in 1:size(sample_t2,2)
        x1_t2 += n_t2[m] / Neff_t2  * sample_t2[:,m] 
        x2_t2 += n_t2[m] / Neff_t2 * sample_t2[:,m] * sample_t2[:,m]'
    end
    r = 0.5
    x1 = r * (x1_t1 + x1_t2) 
    x2 = r * (x2_t1 + x2_t2) 
    
    # Integrated Covariance
    C = copy(x2) - (2 * x1_t1 * x1_t1' + 2 * x1_t2 * x1_t2' + x1_t1 * x1_t2'+ x1_t2 * x1_t1') * (1.0/6)
    
    for i in 1:L
        C[i,i] = (3 - 2 * x1_t2[i]) * (x1_t1[i] + x1_t2[i]) - 2 * x1_t1[i]^2
    end
    #C /= N
    return (x1, x2, C)
end;

#time_scale = "mean" # possible choice should be the "left-edge" and "right-edge". 
function get_Ctot_drift_x1mean(L, data, t_step, time_scale="mean", dg=1.0e-1)
    # treatment of time is naive f(x(t))dx(t) = f(t)*(Dt[t]+Dt[t-1])/2
    # Naive remedy is f(x(t))dx(t) = f(t)*[ (t-t_old)/2 + (t_next-t)/2 ]
    C_tot = zeros(L,L)
    drift = zeros(L)

    x1_ini = zeros(L) 
    x1_fin = zeros(L)
    x1_mean = zeros(L)
    t_accume = 0
    Dt_accume = 0

    length_t_step = length(t_step)
    for i in 1:length_t_step
        t_accume += t_step[i]
        (n_t, sample_t) = get_sample_at_t(data, t_accume); # without interpolate
        (x1,x2,C) = get_x1_x2_C(L, n_t, sample_t); # without interpolate
        Dt = 0
        if(i<length_t_step)
            if(i==1)
                x1_ini = copy(x1)
            end
        end
        if(i==length_t_step)
            x1_fin = copy(x1)
            Dt = t_step[i]
        end

        if(time_scale == "mean")
            if(i<length_t_step)
                Dt = 0.5 * ( t_step[i] + t_step[i+1] )
            end
            if(i==length_t_step)
                Dt = t_step[i]
            end
        end

        if(time_scale == "left-edge")
            Dt = t_step[i]
        end

        if(time_scale == "right-edge")
            if(i<length_t_step)
                Dt = t_step[i+1]
            end
            if(i==length(t_step))
                Dt = t_step[length_t_step]
            end
        end
        Dt_accume += Dt
        drift += dg * Dt * ( ones(L) - 2 * x1)
        C_tot += dg * Dt * C
        x1_mean += x1 * Dt
    end
    x1_mean /= Dt_accume;
    return (C_tot, drift, x1_mean, x1_ini, x1_fin)
end;

selection_coeff(mu, drift, C_tot, x1_ini, x1_fin, r=0.001) = (C_tot + r*I )\(x1_fin - (x1_ini + mu*drift) );

# Example to get the selection coefficient.
#(C_tot, drift, x1_mean, x1_ini, x1_fin) = get_Ctot_drift_x1mean(L, N, data, t_step, "mean", 1.0e-1)
#selection_1 = selction_coeff(mu, drift, C_tot, x1_ini, x1_fin);

function get_integrated_Ctot_drift_x1mean_binary_simple(L, data, time_list_in)
    time_list = unique(time_list_in) 
    C_tot = zeros(L,L)
    drift = zeros(L)
    x1_mean = zeros(L)
    dg = 1.0
    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1_t1, x2_t1) = get_x1_x2(L, n_t, sample_t)
    x1_init = copy(x1_t1)
    x1_mean = copy(x1_t1) * t_old
    Dt_sum = 0 

    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t); # without interpolate
        #@show t, size(sample_t)
        (x1_t2, x2_t2) = get_x1_x2(L, n_t, sample_t)

        (x1_integrated, C_integreated) = get_integrated_corr_mean(t, L, x1_t1, x1_t2, x2_t1, x2_t2)
        x1_t1 = copy(x1_t2); 
        x2_t1 = copy(x2_t2)

         Dt = t - t_old; 
        t_old = t
        
        drift += dg * Dt * (ones(L) - 2 * x1_integrated)
        C_tot += dg * Dt * C_integreated
        x1_mean += Dt * copy(x1_t2) 
        Dt_sum += Dt
    end
    (evl,evt) = eigen(C_tot)
    x1_mean /= (Dt_sum+time_list[1]) 
    x1_fini = copy(x1_t1)
    return (C_tot, drift, x1_mean, x1_init, x1_fini)
end;


function get_naive_Ctot_drift_x1mean(L, data, time_list_in)
    time_list = unique(time_list_in) 
    t_old = time_list[1]
    C_tot = zeros(L,L)
    drift = zeros(L)
    dg = 1.0

    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1, x2) = get_x1_x2(L, n_t, sample_t)
    x1_init = copy(x1)
    x1_mean = copy(x1) * t_old
    Dt_sum = t_old
    x1_fini = zeros(size(x1_init))
    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t)
        (x1, x2) = get_x1_x2(L, n_t, sample_t)
        C = x2 - x1 * x1'
        C_diag = x1 .* (ones(size(x1)) - x1)
        C[diagind(C)] = C_diag
        Dt = t - t_old; t_old = t

        drift += dg * Dt * (ones(L) - 2* x1)
        C_tot += dg * Dt * C
        x1_mean += x1 * Dt
        Dt_sum += Dt

        if(t==time_list[end])
            x1_fini = copy(x1)
        end
    end
    x1_mean /= Dt_sum;
    return (C_tot, drift, x1_mean, x1_init, x1_fini)
end;

