using Pkg
using DelimitedFiles
using Profile    
using Random
using LinearAlgebra
using Plots
using Printf
rng = Random.MersenneTwister(1234);

function eigenp(A)
	(evl,evt) = eigen(A)
	return (abs.(evl), evt)
end

function get_Bezier_Ctot_drift_x1mean_categorical_simple(q, L, data, time_list_in, muMat, pseud_count)
    time_list = unique(time_list_in) 
    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1, x2) = get_x1_x2_psa(q*L, n_t, sample_t, pseud_count)
    #(x1, x2) = get_x1_x2(q, L, n_t, sample_t)
    x1_init = copy(x1)
    
    x2_traject = []; 
    x1_traject = []
    x2_traject = push!(x2_traject, copy(x2))
    x1_traject = push!(x1_traject, copy(x1))
    x1_fini = zeros(size(x1_init))
    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t)
        (x1, x2) = get_x1_x2_psa(q*L, n_t, sample_t, pseud_count)
        #(x1, x2) = get_x1_x2(q, L, n_t, sample_t)
        x2_traject = push!(x2_traject, copy(x2))
        x1_traject = push!(x1_traject, copy(x1))
        if(t==time_list[end])
            x1_fini = copy(x1)
        end
    end

    (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b) = get_Bezier_vec_Mat_for_all_sites_time_categorical(q, L, x1_traject, x2_traject)
    Dt_sum = t_old
    x1_mean = zeros(size(x1_traject[1]))
    B2 = get_integrated_Bezier_2nd()
    
    C_tot1 = zeros(q*L, q*L)
    drift = zeros(q*L)
    for t_id in 1:(length(x1_traject)-1)
        # Integrated single frequencies
        x1_t1 = x1_traject[t_id]; 
        x1_t2 = x1_traject[t_id+1]        
        x1 = 0.25 * ( x1_t1 + point_set_x1_a[:,t_id] + point_set_x1_b[:,t_id] + x1_t2)
            
        # Integrated pairwise frequencies
        x2_t1 = x2_traject[t_id]; 
        x2_t2 = x2_traject[t_id+1]
        x2 = 0.25 * ( x2_t1 + point_set_x2_a[:,:,t_id]+ point_set_x2_b[:,:,t_id] + x2_t2)
        
        # Get integral of frequency products
        B1 = copy(x1_t1)
        B1 = hcat(B1, 3*point_set_x1_a[:,t_id])
        B1 = hcat(B1, 3*point_set_x1_b[:,t_id])
        B1 = hcat(B1, x1_t2)
        
        # Get integrated interpolated covariance matrix
        C_temp1 = B1 * B2 * B1'
        C_temp1 = 0.5 * copy(C_temp1 + C_temp1')
        C1 = x2 - C_temp1
        C_diag1 = copy(x1) - Array([B1[n,:]' * B2 * B1[n,:] for n in 1:size(B1,1)]) 
        C1[diagind(C1)] = C_diag1
        
        t = time_list[t_id+1]
        Dt = t - t_old; 
        t_old = t

        (fluxIn, fluxOut) = get_flux(q, L, x1, muMat) 
        drift += Dt * (fluxIn - fluxOut)
        C_tot1 += Dt * C1
        x1_mean += x1 * Dt
        Dt_sum += Dt
    end
    x1_mean /= Dt_sum

    (evl, evt) = eigen(C_tot1)   # If there are negative eigenvalues, this method will be more unstable. 
    #(evl, evt) = eigenp(C_tot1) # It will be more stable

    return (C_tot1, drift, x1_mean, x1_init, x1_fini, point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;


function get_Bezier_Ctot_drift_x1mean_categorical_simple_poly(
		poly_seq_len, q, L, data, time_list_in, muMat, pseud_count,
		polymorphic_positions, 
		polymorphic_states_set,
		original_index_to_polymorophic_index,
		scaling_psc
		)
    # qL is not necessary to be q*L
    time_list = unique(time_list_in) 
    t_old = time_list[1]
    (n_t, sample_t) = get_sample_at_t(data, t_old)
    (x1, x2) = get_x1_x2_psa(poly_seq_len, n_t, sample_t, scaling_psc, pseud_count)
    #(x1, x2) = get_x1_x2(q, L, n_t, sample_t)
    x1_init = copy(x1)
    x2_traject = []; 
    x1_traject = []
    x2_traject = push!(x2_traject, copy(x2))
    x1_traject = push!(x1_traject, copy(x1))
    x1_fini = zeros(size(x1_init))
    for t in time_list[2:end]
        (n_t, sample_t) = get_sample_at_t(data, t)
	(x1, x2) = get_x1_x2_psa(poly_seq_len, n_t, sample_t, scaling_psc, pseud_count)
        #(x1, x2) = get_x1_x2(q, L, n_t, sample_t)
        x2_traject = push!(x2_traject, copy(x2))
        x1_traject = push!(x1_traject, copy(x1))
        if(t==time_list[end])
            x1_fini = copy(x1)
        end
    end

    (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b) = get_Bezier_vec_Mat_for_all_sites_time_categorical_poly(poly_seq_len, x1_traject, x2_traject)
    # the out put doesn't have q x L structure. 
    Dt_sum = t_old
    x1_mean = zeros(size(x1_traject[1]))
    B2 = get_integrated_Bezier_2nd()
    
    C_tot1 = zeros(poly_seq_len, poly_seq_len)
    drift = zeros(poly_seq_len)
    for t_id in 1:(length(x1_traject)-1)
        # Integrated single frequencies
        x1_t1 = x1_traject[t_id]; 
        x1_t2 = x1_traject[t_id+1]        
        x1 = 0.25 * ( x1_t1 + point_set_x1_a[:,t_id] + point_set_x1_b[:,t_id] + x1_t2)
            
        # Integrated pairwise frequencies
        x2_t1 = x2_traject[t_id]; 
        x2_t2 = x2_traject[t_id+1]
        x2 = 0.25 * ( x2_t1 + point_set_x2_a[:,:,t_id]+ point_set_x2_b[:,:,t_id] + x2_t2)
        
        # Get integral of frequency products
        B1 = copy(x1_t1)
        B1 = hcat(B1, 3*point_set_x1_a[:,t_id])
        B1 = hcat(B1, 3*point_set_x1_b[:,t_id])
        B1 = hcat(B1, x1_t2)
        
        # Get integrated interpolated covariance matrix
        C_temp1 = B1 * B2 * B1'
        C_temp1 = 0.5 * copy(C_temp1 + C_temp1')
        C1 = x2 - C_temp1
        C_diag1 = copy(x1) - Array([B1[n,:]' * B2 * B1[n,:] for n in 1:size(B1,1)]) 
        C1[diagind(C1)] = C_diag1
        
        t = time_list[t_id+1]
        Dt = t - t_old; 
        t_old = t

        (fluxIn, fluxOut) = get_flux_poly(poly_seq_len, q, L, x1, muMat, 
				polymorphic_positions, 
			       	polymorphic_states_set,
			       	original_index_to_polymorophic_index)
	drift += Dt * (fluxIn - fluxOut)
        C_tot1 += Dt * C1
        x1_mean += x1 * Dt
        Dt_sum += Dt
    end

    x1_mean /= Dt_sum

    #(evl, evt) = eigen(C_tot1)   # If there are negative eigenvalues, this method will be more unstable. 
    #(evl, evt) = eigenp(C_tot1) # It will be more stable

    return (C_tot1, drift, x1_mean, x1_init, x1_fini, point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;

##### Suppose each Bezier curve is not independent but is defined as a q dimensional vector.
function get_Bezier_vec_Mat_for_all_sites_time_categorical(q, L, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    #NOTE: index ordering is not q L n, but n L q.  
    point_set_x1 = zeros(L, n_traject, q)
    point_set_x1_a = zeros(q*L, N)
    point_set_x1_b = zeros(q*L, N)
    for i in 1:L
       	for a in 1:q 
            point_set_x1[i,:,a] = [ x1_traject[n][km(i,a,q)] for n in 1:n_traject];
    	end
        (a_vec_i, b_vec_i) = get_a_b_Bezier_categorical_categorical(N, point_set_x1[i,:,:])
        for n in 1:N	
            point_set_x1_a[km.(i,1:q,q),n] = copy(a_vec_i[n,:])
            point_set_x1_b[km.(i,1:q,q),n] = copy(b_vec_i[n,:])
        end    
    end
    
    #-- Get interpolation points for correlations for each time step --#
    point_set_x2 = zeros(L,L, n_traject,q*q)
    point_set_x2_a = zeros(q*L,q*L, N)
    point_set_x2_b = zeros(q*L,q*L, N)
    
    P_tensor = zeros(q*L,q*L, N)
    for n in 2:(N-1)
	    P_tensor[:, :, n] = 2 * ( 2 * x2_traject[n][:,:] + x2_traject[n+1][:,:] ) 
    end
    P_tensor[:,:, 1] = x2_traject[1][:,:] + 2 * x2_traject[2][:,:]  
    P_tensor[:,:, N] = 8 * x2_traject[N][:,:] + x2_traject[N+1][:,:]  

    Bez_mat = get_BezMat(N)
    inv_Bez = inv(Bez_mat) 

    for n in 1:N
	    for m in 1:N
		    point_set_x2_a[:,:,n] += inv_Bez[n, m] * P_tensor[:, :, m] 
	    end
    end
    for n in 1:(N-1)
	    point_set_x2_b[:,:,n] = 2 * x2_traject[n+1][:,:] - point_set_x2_a[:,:,n+1]
    end
    point_set_x2_b[:,:,N] = 0.5 * (point_set_x2_a[:,:,N] + x2_traject[N+1][:,:]) 

    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;

function selction_coeff_cholfact(mu, reg, drift, C_tot, x1_init, x1_fini)
    return s = (C_tot+reg*I) \ (x1_fini-x1_init-mu*drift)
end;

##### Suppose each Bezier curve is not independent but is defined as a q dimensional vector.
function get_Bezier_vec_Mat_for_all_sites_time_categorical(q, L, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    #NOTE: index ordering is not q L n, but n L q.  
    point_set_x1 = zeros(L, n_traject, q)
    point_set_x1_a = zeros(q*L, N)
    point_set_x1_b = zeros(q*L, N)
    for i in 1:L
       	for a in 1:q 
            point_set_x1[i,:,a] = [ x1_traject[n][km(i,a,q)] for n in 1:n_traject];
    	end
        (a_vec_i, b_vec_i) = get_a_b_Bezier_categorical_categorical(N, point_set_x1[i,:,:])
        for n in 1:N	
            point_set_x1_a[km.(i,1:q,q),n] = copy(a_vec_i[n,:])
            point_set_x1_b[km.(i,1:q,q),n] = copy(b_vec_i[n,:])
        end    
    end
    
    #-- Get interpolation points for correlations for each time step --#
    point_set_x2 = zeros(L,L, n_traject,q*q)
    point_set_x2_a = zeros(q*L,q*L, N)
    point_set_x2_b = zeros(q*L,q*L, N)
    
    P_tensor = zeros(q*L,q*L, N)
    for n in 2:(N-1)
	    P_tensor[:, :, n] = 2 * ( 2 * x2_traject[n][:,:] + x2_traject[n+1][:,:] ) 
    end
    P_tensor[:,:, 1] = x2_traject[1][:,:] + 2 * x2_traject[2][:,:]  
    P_tensor[:,:, N] = 8 * x2_traject[N][:,:] + x2_traject[N+1][:,:]  

    Bez_mat = get_BezMat(N)
    inv_Bez = inv(Bez_mat) 

    for n in 1:N
	    for m in 1:N
		    point_set_x2_a[:,:,n] += inv_Bez[n, m] * P_tensor[:, :, m] 
	    end
    end
    for n in 1:(N-1)
	    point_set_x2_b[:,:,n] = 2 * x2_traject[n+1][:,:] - point_set_x2_a[:,:,n+1]
    end
    point_set_x2_b[:,:,N] = 0.5 * (point_set_x2_a[:,:,N] + x2_traject[N+1][:,:]) 

    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;

##### Suppose each Bezier curve is not independent but is defined as a q dimensional vector.
function get_Bezier_vec_Mat_for_all_sites_time_categorical_poly(qL, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    #NOTE: index ordering is not q L n, but n L q.  
    point_set_x1 = zeros(qL, n_traject)
    point_set_x1_a = zeros(qL, N)
    point_set_x1_b = zeros(qL, N)
    for i in 1:qL
        point_set_x1[i,:] = [ x1_traject[n][i] for n in 1:n_traject];
        (a_vec_i, b_vec_i) = get_a_b_Bezier_categorical_categorical(N, point_set_x1[i,:])
        point_set_x1_a[i,:] = copy(a_vec_i)
        point_set_x1_b[i,:] = copy(b_vec_i)
    end
    
    #-- Get interpolation points for correlations for each time step --#
    point_set_x2 = zeros(qL,qL, n_traject)
    point_set_x2_a = zeros(qL,qL, N)
    point_set_x2_b = zeros(qL,qL, N)
    
    P_tensor = zeros(qL,qL, N)
    for n in 2:(N-1)
	    P_tensor[:, :, n] = 2 * ( 2 * x2_traject[n][:,:] + x2_traject[n+1][:,:] ) 
    end
    P_tensor[:,:, 1] = x2_traject[1][:,:] + 2 * x2_traject[2][:,:]  
    P_tensor[:,:, N] = 8 * x2_traject[N][:,:] + x2_traject[N+1][:,:]  

    Bez_mat = get_BezMat(N)
    inv_Bez = inv(Bez_mat) 

    for n in 1:N
	    for m in 1:N
		    point_set_x2_a[:,:,n] += inv_Bez[n, m] * P_tensor[:, :, m] 
	    end
    end
    for n in 1:(N-1)
	    point_set_x2_b[:,:,n] = 2 * x2_traject[n+1][:,:] - point_set_x2_a[:,:,n+1]
    end
    point_set_x2_b[:,:,N] = 0.5 * (point_set_x2_a[:,:,N] + x2_traject[N+1][:,:]) 
    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;

# a_vec, b_vec should be N x D dimension 
function get_a_b_Bezier_categorical_categorical(N, p_set)
    Bez_mat = get_BezMat(N)
    Bez_vec = get_BezVec_categorical(N, p_set);
    
    #Ax=b
    #x = cholesky(A) \ b # if A is Hermetian
    a_vec = Bez_mat \ Bez_vec;
    b_vec = get_a2b_Bez_categorical(p_set, a_vec);
    return (a_vec, b_vec)
end;

"""
This is 3rd order Bezier
 x is N+1 x D (=q or q x q) matrix
 a is N x D (=q or q x q) matrix
"""
function get_a2b_Bez_categorical(x,a)
	(N,D) = size(a)
	b = zeros(N,D)
	for n in 1:(N-1)
	    b[n,:] = 2*x[n+1,:] - a[n+1,:]
	end
	b[N,:] = 0.5*(a[N,:] + x[N+1,:])
	return b
end;

function get_BezMat(N)
    Bez_mat = zeros(N,N);
    for n in 2:(N-1)
        if(n>1 && n<N)
            Bez_mat[n,n] = 4; 
	    Bez_mat[n,n+1] = 1; 
	    Bez_mat[n,n-1] = 1;
        end
    end
    Bez_mat[1,1] = 2; 
    Bez_mat[1,2] = 1;
    Bez_mat[N,N-1] = 2; 
    Bez_mat[N,N] = 7;
    return Bez_mat
end;

# x is N x D(=q, =q^2) matrix
function get_BezVec_categorical(N, x)
	D = size(x,2)
	Bez_vec = zeros(N,D);
	for n in 2:(N-1)
	    Bez_vec[n,:] = 2*(2*x[n,:] + x[n+1,:])
	end
	Bez_vec[1,:] = x[1,:] + 2*x[2,:]
	Bez_vec[N,:] = 8*x[N,:] + x[N+1,:]
	return Bez_vec
end;

function get_integrated_Bezier_2nd()
    a = 1.0/7
    b = 1.0/42
    c = 1.0/105
    d = 1.0/140
    return B2 = Matrix([[a, b, c, d] [b, c, d, c] [c, d, c, b] [d, c, b, a]])
end;
