# suppose a_vec, b_vec are N x q dimmensions
function filter_controle_points(L, a_vec, b_vec, p_set)
    for n in 1:size(a_vec,1)
	if(abs(p_set[n]-p_set[n+1])<1e-2)
	    p_mean = 0.5 * (p_set[n]+p_set[n+1]) 
            a_vec[n] = p_mean
            b_vec[n] = p_mean
        end
        if(a_vec[n]<0)
            a_vec[n] = 0      
        end
        if(b_vec[n]<0)
            b_vec[n] = 0      
        end
    end
    return (a_vec, b_vec)
end;

function get_Bezier_vec_Mat_for_all_sites_time_binary_filter(L, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject, 1)
    N = n_traject - 1
    
    point_set_x1_a = zeros(L, N)
    point_set_x1_b = zeros(L, N)
    point_set_x2_a = zeros(L, L, N)
    point_set_x2_b = zeros(L, L, N)
     
    P2_tensor = zeros(L,L, N)
    P1_tensor = zeros(L, N)	
    for n in 2:(N-1)
	P2_tensor[:, :, n] = 2 * ( 2 * x2_traject[n][:,:] + x2_traject[n+1][:,:] ) 
	P1_tensor[:, n] = 2 * ( 2 * x1_traject[n][:] + x1_traject[n+1][:] ) 
    end
    P2_tensor[:,:, 1] = x2_traject[1][:,:] + 2 * x2_traject[2][:,:]  
    P2_tensor[:,:, N] = 8 * x2_traject[N][:,:] + x2_traject[N+1][:,:]  
    P1_tensor[:, 1] = x1_traject[1][:] + 2 * x1_traject[2][:]  
    P1_tensor[:, N] = 8 * x1_traject[N][:] + x1_traject[N+1][:]  

    Bez_mat = get_BezMat(N)
    inv_Bez = inv(Bez_mat) 

    for n in 1:N
	    for m in 1:N
		    point_set_x2_a[:,:,n] += inv_Bez[n, m] * P2_tensor[:, :, m] 
		    point_set_x1_a[:,n] += inv_Bez[n, m] * P1_tensor[:, m] 
	    end
    end
    for n in 1:(N-1)
	    point_set_x2_b[:,:,n] = 2 * x2_traject[n+1][:,:] - point_set_x2_a[:,:,n+1]
	    point_set_x1_b[:,n] = 2 * x1_traject[n+1][:] - point_set_x1_a[:,n+1]
    end
    point_set_x2_b[:,:,N] = 0.5 * (point_set_x2_a[:,:,N] + x2_traject[N+1][:,:]) 
    point_set_x1_b[:,N] = 0.5 * (point_set_x1_a[:,N] + x1_traject[N+1][:]) 
   
    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;


function get_Bezier_Ctot_drift_x1mean_binary_simple(L, data, time_list_in)
    time_list = unique(time_list_in) 
    dg = 1.0

    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1, x2) = get_x1_x2(L, n_t, sample_t)
    x1_init = copy(x1)
    
    x2_traject = []; 
    x1_traject = []
    x2_traject = push!(x2_traject, copy(x2))
    x1_traject = push!(x1_traject, copy(x1))
    x1_fini = zeros(size(x1_init))
    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t)
        (x1, x2) = get_x1_x2(L, n_t, sample_t)
        x2_traject = push!(x2_traject, copy(x2))
        x1_traject = push!(x1_traject, copy(x1))
        if(t==time_list[end])
            x1_fini = copy(x1)
        end
    end

    (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b) = get_Bezier_vec_Mat_for_all_sites_time_binary_filter(L, x1_traject, x2_traject)

    Dt_sum = t_old
    x1_mean = zeros(size(x1_traject[1]))
    B2 = get_integrated_Bezier_2nd()
    
    C_tot1 = zeros(L, L); C_tot2 = zeros(L, L)
    drift = zeros(L)
    for t_id in 1:(length(x1_traject)-1)
        # Integrated single freq.
        x1_t1 = x1_traject[t_id]; 
        x1_t2 = x1_traject[t_id+1]
        x1 = 0.25 * ( x1_t1 + point_set_x1_a[:,t_id] + point_set_x1_b[:,t_id] + x1_t2)
            
        # Integrated single freq.
        x2_t1 = x2_traject[t_id]; 
        x2_t2 = x2_traject[t_id+1]
    
        x2 = 0.25 * ( x2_t1 + point_set_x2_a[:,:,t_id]+ point_set_x2_b[:,:,t_id] + x2_t2)
        
        B1 = copy(x1_t1)
        B1 = hcat(B1, 3*point_set_x1_a[:,t_id])
        B1 = hcat(B1, 3*point_set_x1_b[:,t_id])
        B1 = hcat(B1, x1_t2)
        
        C_temp1 = B1 * B2 * B1'
        C_temp1 = 0.5 * copy(C_temp1 + C_temp1')
        C_temp2 = x1 * x1'
        
        C1 = x2 - C_temp1
        C_diag1 = copy(x1) - Array([B1[n,:]' * B2 * B1[n,:] for n in 1:size(B1,1)]) 
        C1[diagind(C1)] = C_diag1
        
        t = time_list[t_id+1]
        Dt = t - t_old; 
        t_old = t
        
        drift += dg * Dt * (ones(L) - 2 * x1)
        C_tot1 += dg * Dt * C1
        
        x1_mean += x1 * Dt
        Dt_sum += Dt
    end
    
    (evl1,evt1) = eigen(C_tot1)
    #@printf("min(λtot1):%.4f, sum(λtot1):%.4f, min(λtot2):%.4f, sum(λtot2):%.4f, sumDt:%d \n\n", minimum(evl1), sum(evl1), minimum(evl2), sum(evl2), Dt_sum)
    
    
    x1_mean /= Dt_sum
    return (C_tot1, drift, x1_mean, x1_init, x1_fini, point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;
