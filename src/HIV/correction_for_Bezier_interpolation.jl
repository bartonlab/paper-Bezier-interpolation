function get_polymorophic_index_set(q, L, polymorphic_sites)
    polymorphic_positions = []
    polymorphic_states_set = []
    original_index_to_polymorophic_index = []
    count_poly_index = 1
    for i in 1:L
        poly_set_i = polymorphic_sites[km.(i,1:q,q) .+ 3]
        if(count(poly_set_i)>1) # count(poly_set_i) will be 0 or >1, (1 will not comes here. )
            push!(polymorphic_positions, i)
            poly_sets_set_i = []
            for a in 1:q
                if(poly_set_i[a])
                    push!(poly_sets_set_i,a)
                end

                #------- make index mapping --------#
                flag_poly = polymorphic_sites[km(i,a,q) .+ 3]
                if(flag_poly)
                    push!(original_index_to_polymorophic_index, count_poly_index)
                    count_poly_index += 1
                end
                if(!flag_poly)
                    push!(original_index_to_polymorophic_index, -1)
                end
                #------------------------------------#
            end
            push!(polymorphic_states_set, copy(poly_sets_set_i)) 
        end
        # Purely monomorphic sites 
        if(count(poly_set_i)==0)
            for a in 1:q
                push!(original_index_to_polymorophic_index, -1)
            end
        end
    end;
    return (
        polymorphic_positions, 
        polymorphic_states_set, 
        original_index_to_polymorophic_index)
end;

function get_polymorphic_scaling_for_pseudo_count(polymorphic_states_set)
    # alpha_set is polymorphic scaling for pseudo count.
    alpha_set = []
    for x in polymorphic_states_set
        n_in = length(x)
        scale_in = 1.0/n_in
        for _ in x
            push!(alpha_set, scale_in)
        end
    end
    return alpha_set
end;


function get_selection(reg, num, C_tot)
    s = (C_tot + reg*I) \ num
end;

function get_selection_SL(reg, num, C_tot)
    s = num ./ (C_tot[diagind(C_tot)] .+ reg)
end;

function write_MPL_SL(est_MPL, est_SL, fnameout_MPL, fnameout_SL)
    fout_MPL = open(fnameout_MPL, "w")
    fout_SL = open(fnameout_SL, "w")
    for x in est_MPL
        println(fout_MPL, x)
    end
    for x in est_SL
        println(fout_SL, x)
    end
    close(fout_MPL); close(fout_SL)
    
end;


# If gaps becoems larger than dt_max, include an imaginary point inbeteween the points. 

function get_Bezier_Ctot_drift_x1mean_categorical_simple_poly_virtual_insertion(
		poly_seq_len, q, L, data, time_list_in, muMat, pseud_count,
		polymorphic_positions, 
		polymorphic_states_set,
		original_index_to_polymorophic_index,
		scaling_psc, 
        time_interval_threshold
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
    (time_list, x1_traject, x2_traject) = get_inserted_time_trajectory_set(time_list, x1_traject, x2_traject, time_interval_threshold)

    
    (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b) = get_Bezier_vec_Mat_for_all_sites_time_categorical_poly(poly_seq_len, x1_traject, x2_traject)
    # the out put doesn't have q x L structure. 
    Dt_sum = t_old
    x1_mean = zeros(size(x1_traject[1]))
    B2 = get_integrated_Bezier_2nd()

    C_tot1 = zeros(poly_seq_len, poly_seq_len)
    C_tot1_positive = zeros(poly_seq_len, poly_seq_len)
    C_tot1_negative = zeros(poly_seq_len, poly_seq_len)
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

        C_tot1_positive += Dt * x2
        C_tot1_negative += Dt * C_temp1
        
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

    return (C_tot1, drift, x1_mean, x1_init, x1_fini, 
            point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b, 
            x1_traject, x2_traject, time_list, 
            C_tot1_positive, C_tot1_negative)
end;

function get_inserted_time_trajectory_set(time_list, x1_traject, x2_traject, time_interval_threshold=50)
    time_list_inserted = []
    x1_traject_inserted = []
    x2_traject_inserted = []
    len_time_list = length(time_list)
    for i_t in 1:(len_time_list-1)
        push!(time_list_inserted, time_list[i_t])
        push!(x1_traject_inserted, copy(x1_traject[i_t]))
        push!(x2_traject_inserted, copy(x2_traject[i_t]))

        t_interval = time_list[i_t+1] - time_list[i_t]
        if(t_interval>time_interval_threshold)
            intermediate_time = 0.5*(time_list[i_t+1] + time_list[i_t])
            push!(time_list_inserted, intermediate_time)

            #also make interpolation for the frequncies. 
            x1_traject_intermediate = copy(0.5*(x1_traject[i_t+1] + x1_traject[i_t]))
            x2_traject_intermediate = copy(0.5*(x2_traject[i_t+1] + x2_traject[i_t]))
            push!(x1_traject_inserted, copy(x1_traject_intermediate))
            push!(x2_traject_inserted, copy(x2_traject_intermediate))
        end
    end;
    push!(time_list_inserted, time_list[end])
    push!(x1_traject_inserted, copy(x1_traject[end]))
    push!(x2_traject_inserted, copy(x2_traject[end]))
    
    return (time_list_inserted, x1_traject_inserted, x2_traject_inserted)
end;


# This function compute the scalar controle points 
function get_Bezier_vec_Mat_for_all_sites_time_poly_element_wise(x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    #NOTE: index ordering is not q L n, but n L q.  

    point_set_x1 = [ x1_traject[n] for n in 1:n_traject];
    (a_vec_i, b_vec_i) = get_a_b_Bezier_categorical_categorical(N, point_set_x1)
    point_set_x1_a = copy(a_vec_i)
    point_set_x1_b = copy(b_vec_i)

    
    #-- Get interpolation points for correlations for each time step --#
    point_set_x2 = zeros(n_traject)
    point_set_x2_a = zeros(N)
    point_set_x2_b = zeros(N)
    
    P_tensor = zeros(N)
    for n in 2:(N-1)
	    P_tensor[n] = 2 * ( 2 * x2_traject[n] + x2_traject[n+1]) 
    end
    P_tensor[1] = x2_traject[1] + 2 * x2_traject[2]  
    P_tensor[N] = 8 * x2_traject[N] + x2_traject[N+1]

    Bez_mat = get_BezMat(N)
    inv_Bez = inv(Bez_mat) 

    for n in 1:N
	    for m in 1:N
		    point_set_x2_a[n] += inv_Bez[n, m] * P_tensor[m] 
	    end
    end
    for n in 1:(N-1)
	    point_set_x2_b[n] = 2 * x2_traject[n+1] - point_set_x2_a[n+1]
    end
    point_set_x2_b[N] = 0.5 * (point_set_x2_a[N] + x2_traject[N+1]) 
    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;


function correction_of_covariance_diagonal_with_midpoint(
        C_tot,
        poly_seq_len, 
        x1_traject, 
        x2_traject,
        time_list,
        diff_freq_thresh=0.6
    )
    
    unique_collecting_date = sort(unique(time_list))
    C_diag_with_mid_point = zeros(poly_seq_len)
    #i_poly = 1 # polymorphic index  
    for i_poly in 1:poly_seq_len
        C_diag_with_mid_point[i_poly] = C_tot[i_poly, i_poly]
        
        (flag_need_to_insert, x1_traject_i, x2_traject_ii) = 
                check_need_to_inclue_middpoint_element_wise(x1_traject, x2_traject, i_poly, diff_freq_thresh)

        if(flag_need_to_insert)
            # the following process should be executed if flag_need_to_insert == true

            (x1_traject_i_with_midpts, x2_traject_ii_with_midpts, unique_collecting_date_with_midpts) = 
                    get_midpoint_freq_element_wise_include_neighbor(x1_traject_i, x2_traject_ii, unique_collecting_date, 
                                                diff_freq_thresh)
            
            # Get controle points 
            # Here, effectively the operation is same as the case where q and L are 1 and 1, respectively.
            #@show x1_traject_i_with_midpts
            (point_set_x1_a_B_i, point_set_x1_b_B_i, point_set_x2_a_B_i, point_set_x2_b_B_i) = 
            get_Bezier_vec_Mat_for_all_sites_time_poly_element_wise(x1_traject_i_with_midpts, x2_traject_ii_with_midpts);

            #Effectively this cell changes the only diagonal parts.
            # the out put doesn't have q x L structure. 
            t_old = unique_collecting_date_with_midpts[1]
            Dt_sum = t_old
            x1_mean_elem = 0.0
            B2 = get_integrated_Bezier_2nd()
            
            Ctot_elem = 0.0 
            #drift_elem = 0.0
            for t_id in 1:(length(x1_traject_i_with_midpts)-1)
                # Integrated single frequencies
                x1_t1_elem = x1_traject_i_with_midpts[t_id]; 
                x1_t2_elem = x1_traject_i_with_midpts[t_id+1]        
                x1_elem = 0.25 * ( x1_t1_elem + point_set_x1_a_B_i[t_id] + point_set_x1_b_B_i[t_id] + x1_t2_elem)

                # Integrated pairwise frequencies
                x2_t1_elem = x2_traject_ii_with_midpts[t_id]; 
                x2_t2_elem = x2_traject_ii_with_midpts[t_id+1]
                x2_elem = 0.25 * ( x2_t1_elem + point_set_x2_a_B_i[t_id]+ point_set_x2_b_B_i[t_id] + x2_t2_elem)

                # Get integral of frequency products
                B1=[]
                push!(B1, x1_t1_elem)
                push!(B1, 3*point_set_x1_a_B_i[t_id])
                push!(B1, 3*point_set_x1_b_B_i[t_id])
                push!(B1, x1_t2_elem)
                B1 = vec(B1)

                # Get integrated interpolated covariance matrix
                C_temp1_elem = B1' * B2 * B1
                C1_elem = x2_elem - C_temp1_elem

                t = unique_collecting_date_with_midpts[t_id+1]
                Dt = t - t_old; 
                t_old = t
                Ctot_elem += Dt * C1_elem
                x1_mean_elem += x1_elem * Dt
                Dt_sum += Dt
            end
            C_diag_with_mid_point[i_poly] = Ctot_elem
        end;
    end
    
    return C_diag_with_mid_point
end;


function correction_of_covariance_diagonal_with_midpoint_on_off_diagonal(
        C_tot,
        C_tot_positive, 
        C_tot_negative,
        poly_seq_len, 
        x1_traject, 
        x2_traject,
        time_list,
        diff_freq_thresh=0.6
    )
    
    unique_collecting_date = sort(unique(time_list))
    C_with_mid_point = copy(C_tot)
    
    C_with_mid_point_positive = zeros(poly_seq_len, poly_seq_len)
    C_with_mid_point_negative = zeros(poly_seq_len, poly_seq_len)
    Matrix_modified_element = zeros(poly_seq_len, poly_seq_len) # give 0 non modified, give 1 modified
    for i_poly in 1:poly_seq_len
        for j_poly in (i_poly+1):poly_seq_len
            #@show i_poly, j_poly
            (flag_need_to_insert, x1_traject_i, x1_traject_j, x2_traject_ij) = 
                    check_need_to_inclue_middpoint_element_wise_on_off_diag(
                        x1_traject, x2_traject, i_poly, j_poly, diff_freq_thresh)

            if(flag_need_to_insert)
            #if(true)
                # the following process should be executed if flag_need_to_insert == true
                (x1_traject_i_with_midpts, x1_traject_j_with_midpts, x2_traject_ij_with_midpts, 
                        unique_collecting_date_with_midpts) = 
                        get_midpoint_freq_element_wise_include_neighbor_on_off_diag(
                            x1_traject_i, x1_traject_j, x2_traject_ij, unique_collecting_date, diff_freq_thresh)
            
                # Get controle points 
                (point_set_x1_a_B_i, point_set_x1_b_B_i, 
                    point_set_x1_a_B_j, point_set_x1_b_B_j, 
                    point_set_x2_a_B_ij, point_set_x2_b_B_ij) = 
                        get_Bezier_vec_Mat_for_all_sites_time_poly_element_wise_on_off_diag(
                            x1_traject_i_with_midpts, x1_traject_j_with_midpts, x2_traject_ij_with_midpts);

                #Effectively this cell changes the only diagonal parts.
                t_old = unique_collecting_date_with_midpts[1]
                Dt_sum = t_old
                B2 = get_integrated_Bezier_2nd()

                Ctot_elem_ij_positive = 0.0 # this is for pairwise-frequency
                Ctot_elem_ij_negative = 0.0 # this is for rank 1 matrix 
                for t_id in 1:(length(x1_traject_i_with_midpts)-1)
                    # Integrated single frequencies
                    x1_t1_elem_i = x1_traject_i_with_midpts[t_id]; 
                    x1_t2_elem_i = x1_traject_i_with_midpts[t_id+1]                    
                    x1_elem_i = 0.25 * ( x1_t1_elem_i + point_set_x1_a_B_i[t_id] + point_set_x1_b_B_i[t_id] + x1_t2_elem_i)

                    x1_t1_elem_j = x1_traject_j_with_midpts[t_id]; 
                    x1_t2_elem_j = x1_traject_j_with_midpts[t_id+1]
                    x1_elem_j = 0.25 * ( x1_t1_elem_j + point_set_x1_a_B_j[t_id] + point_set_x1_b_B_j[t_id] + x1_t2_elem_j)
                   
                    # Integrated pairwise frequencies
                    x2_t1_elem_ij = x2_traject_ij_with_midpts[t_id]; 
                    x2_t2_elem_ij = x2_traject_ij_with_midpts[t_id+1]
                    x2_elem_ij = 0.25 * ( x2_t1_elem_ij + point_set_x2_a_B_ij[t_id]+ point_set_x2_b_B_ij[t_id] + x2_t2_elem_ij)

                    # Get integral of frequency products
                    B1_i=[]
                    push!(B1_i, x1_t1_elem_i)
                    push!(B1_i, 3*point_set_x1_a_B_i[t_id])
                    push!(B1_i, 3*point_set_x1_b_B_i[t_id])
                    push!(B1_i, x1_t2_elem_i)
                    B1_i = vec(B1_i)
                    
                    B1_j=[]
                    push!(B1_j, x1_t1_elem_j)
                    push!(B1_j, 3*point_set_x1_a_B_j[t_id])
                    push!(B1_j, 3*point_set_x1_b_B_j[t_id])
                    push!(B1_j, x1_t2_elem_j)
                    B1_j = vec(B1_j)
                    

                    # Get integrated interpolated covariance matrix
                    C_temp1_elem_ij = B1_i' * B2 * B1_j
                    C1_elem_ij = x2_elem_ij - C_temp1_elem_ij

                    t = unique_collecting_date_with_midpts[t_id+1]
                    Dt = t - t_old; 
                    t_old = t

                    Ctot_elem_ij_positive += Dt * x2_elem_ij
                    Ctot_elem_ij_negative += Dt * C_temp1_elem_ij
                    
                    #Ctot_elem_ij += Dt * C1_elem_ij
                    Dt_sum += Dt
                end
                C_with_mid_point[i_poly, j_poly] = Ctot_elem_ij_positive - Ctot_elem_ij_negative
                C_with_mid_point[j_poly, i_poly] = Ctot_elem_ij_positive - Ctot_elem_ij_negative
                #@show Ctot_elem_ij_positive
                C_with_mid_point_positive[j_poly, i_poly] = Ctot_elem_ij_positive
                C_with_mid_point_positive[i_poly, j_poly] = Ctot_elem_ij_positive
                C_with_mid_point_negative[j_poly, i_poly] = Ctot_elem_ij_negative
                C_with_mid_point_negative[i_poly, j_poly] = Ctot_elem_ij_negative
                Matrix_modified_element[i_poly, j_poly] = 1
                Matrix_modified_element[j_poly, i_poly] = 1
            end
            
        end;
    end
    
    return (C_with_mid_point, C_with_mid_point_positive, C_with_mid_point_negative, Matrix_modified_element)
end;


# suppose the i and j are polymorophic indices.
# Note this function need to retrun flag_inserted=true even if encount only one time interval. 
function check_need_to_inclue_middpoint_element_wise_on_off_diag(
        x1_traject, x2_traject, i, j, diff_freq_thresh=0.6)
    flag_inserted = false
    #if frequency change more than 0.6 % then make a frag
    idx_t_max = size(x1_traject, 1)
    x1_traject_i = []
    x1_traject_j = []
    x2_traject_ij = []
    
    push!(x1_traject_i, copy(x1_traject[1][i]))
    push!(x1_traject_j, copy(x1_traject[1][j]))
    push!(x2_traject_ij, copy(x2_traject[1][i, j]))
    for i_t in 2:idx_t_max
        push!(x1_traject_i, copy(x1_traject[i_t][i]))
        push!(x1_traject_j, copy(x1_traject[i_t][j]))
        push!(x2_traject_ij, copy(x2_traject[i_t][i,j]))
        for a in 1:q
            diff_freq_i = abs(x1_traject[i_t][i] - x1_traject[i_t-1][i])
            diff_freq_j = abs(x1_traject[i_t][j] - x1_traject[i_t-1][j])
            diff_freq_ij = abs(x2_traject[i_t][i,j] - x2_traject[i_t-1][i,j])
            
            if(diff_freq_i>=diff_freq_thresh || diff_freq_j>=diff_freq_thresh || diff_freq_ij>=diff_freq_thresh)
                flag_inserted = true
                break;
            end
        end
    end
    return (flag_inserted, x1_traject_i, x1_traject_j, x2_traject_ij)
end;



function get_midpoint_freq_element_wise_include_neighbor_on_off_diag(
        x1_traject_i, x1_traject_j, x2_traject_ij, unique_collecting_date, diff_freq_thresh = 0.6)
    flag_inserted = false
    idx_t_max = size(x1_traject_i, 1)
    midpoint_vec = []
    x1_traject_i_with_midpts = []
    x1_traject_j_with_midpts = []
    x2_traject_ij_with_midpts = []
    unique_collecting_date_with_midpts = []
    
    push!(x1_traject_i_with_midpts, copy(x1_traject_i[1]))
    push!(x1_traject_j_with_midpts, copy(x1_traject_j[1]))
    push!(x2_traject_ij_with_midpts, copy(x2_traject_ij[1]))
    push!(unique_collecting_date_with_midpts, unique_collecting_date[1])
    
    mixing_rate = 0.999 # this is rato for the actually observed frequency
    
    for i_t in 2:idx_t_max
        diff_freq_i = abs(x1_traject_i[i_t] - x1_traject_i[i_t-1])
        diff_freq_j = abs(x1_traject_j[i_t] - x1_traject_j[i_t-1])
        if(diff_freq_i>=diff_freq_thresh || diff_freq_j>=diff_freq_thresh)
            # these are not vector but scalar variable.
            midpoint_i = 0.5*(x1_traject_i[i_t] + x1_traject_i[i_t-1])
            midpoint_j = 0.5*(x1_traject_j[i_t] + x1_traject_j[i_t-1])
            midpoint_ij = 0.5*(x2_traject_ij[i_t] + x2_traject_ij[i_t-1])
            midpoint_time = 0.5*(unique_collecting_date[i_t] + unique_collecting_date[i_t-1])
            push!(x1_traject_i_with_midpts, copy(midpoint_i))
            push!(x1_traject_j_with_midpts, copy(midpoint_j))
            push!(x2_traject_ij_with_midpts, copy(midpoint_ij))
            push!(unique_collecting_date_with_midpts, copy(midpoint_time))
            flag_inserted = true
            
            #insert midpoint before this interval if it defined more than 2-3th points
            if(i_t>=3)
                midpoint_i = mixing_rate*x1_traject_i[i_t-1] + (1-mixing_rate)*x1_traject_i[i_t-2]
                midpoint_j = mixing_rate*x1_traject_j[i_t-1] + (1-mixing_rate)*x1_traject_j[i_t-2]
                midpoint_ij = mixing_rate*x2_traject_ij[i_t-1] + (1-mixing_rate)*x2_traject_ij[i_t-2]
                midpoint_time = mixing_rate*unique_collecting_date[i_t-1] + (1-mixing_rate)*unique_collecting_date[i_t-2]
                push!(x1_traject_i_with_midpts, copy(midpoint_i))
                push!(x1_traject_j_with_midpts, copy(midpoint_j))
                push!(x2_traject_ij_with_midpts, copy(midpoint_ij))
                push!(unique_collecting_date_with_midpts, copy(midpoint_time))                
            end
            
            #insert midpoint after this interval if neighboring interval i_t+1 <= idx_t_max
            if(i_t+1 <= idx_t_max)
                midpoint_i = (1-mixing_rate)*x1_traject_i[i_t+1] + mixing_rate*x1_traject_i[i_t]
                midpoint_j = (1-mixing_rate)*x1_traject_j[i_t+1] + mixing_rate*x1_traject_j[i_t]
                midpoint_ij = (1-mixing_rate)*x2_traject_ij[i_t+1] + mixing_rate*x2_traject_ij[i_t]
                midpoint_time = (1-mixing_rate)*unique_collecting_date[i_t+1] + mixing_rate*unique_collecting_date[i_t]
                push!(x1_traject_i_with_midpts, copy(midpoint_i))
                push!(x1_traject_j_with_midpts, copy(midpoint_j))
                push!(x2_traject_ij_with_midpts, copy(midpoint_ij))
                push!(unique_collecting_date_with_midpts, copy(midpoint_time))                
            end
        end        
        push!(x1_traject_i_with_midpts, copy(x1_traject_i[i_t]))
        push!(x1_traject_j_with_midpts, copy(x1_traject_j[i_t]))
        push!(x2_traject_ij_with_midpts, copy(x2_traject_ij[i_t]))
        push!(unique_collecting_date_with_midpts, copy(unique_collecting_date[i_t]))
    end
    
    index_sorting = sortperm(unique_collecting_date_with_midpts)
    x1_traject_i_with_midpts = copy(x1_traject_i_with_midpts[index_sorting])
    x1_traject_j_with_midpts = copy(x1_traject_j_with_midpts[index_sorting])
    x2_traject_ij_with_midpts = copy(x2_traject_ij_with_midpts[index_sorting])
    unique_collecting_date_with_midpts = copy(unique_collecting_date_with_midpts[index_sorting])
    
    return (x1_traject_i_with_midpts, x1_traject_j_with_midpts, 
        x2_traject_ij_with_midpts, unique_collecting_date_with_midpts)
end;



# This function compute the scalar controle points 
function get_Bezier_vec_Mat_for_all_sites_time_poly_element_wise_on_off_diag(
        x1_traject_i, x1_traject_j, x2_traject_ij)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject_i,1)
    N = n_traject-1
    #NOTE: index ordering is not q L n, but n L q.  

    point_set_x1_i = [ x1_traject_i[n] for n in 1:n_traject];
    (a_vec_i_i, b_vec_i_i) = get_a_b_Bezier_categorical_categorical(N, point_set_x1_i)
    point_set_x1_a_i = copy(a_vec_i_i)
    point_set_x1_b_i = copy(b_vec_i_i)

    point_set_x1_j = [ x1_traject_j[n] for n in 1:n_traject];
    (a_vec_i_j, b_vec_i_j) = get_a_b_Bezier_categorical_categorical(N, point_set_x1_j)
    point_set_x1_a_j = copy(a_vec_i_j)
    point_set_x1_b_j = copy(b_vec_i_j)
    
    #-- Get interpolation points for correlations for each time step --#
    point_set_x2_ij = zeros(n_traject)
    point_set_x2_a_ij = zeros(N)
    point_set_x2_b_ij = zeros(N)
    
    P_tensor_ij = zeros(N)
    for n in 2:(N-1)
	    P_tensor_ij[n] = 2 * ( 2 * x2_traject_ij[n] + x2_traject_ij[n+1]) 
    end
    P_tensor_ij[1] = x2_traject_ij[1] + 2 * x2_traject_ij[2]  
    P_tensor_ij[N] = 8 * x2_traject_ij[N] + x2_traject_ij[N+1]

    Bez_mat = get_BezMat(N)
    inv_Bez = inv(Bez_mat) 

    for n in 1:N
	    for m in 1:N
		    point_set_x2_a_ij[n] += inv_Bez[n, m] * P_tensor_ij[m] 
	    end
    end
    for n in 1:(N-1)
	    point_set_x2_b_ij[n] = 2 * x2_traject_ij[n+1] - point_set_x2_a_ij[n+1]
    end
    point_set_x2_b_ij[N] = 0.5 * (point_set_x2_a_ij[N] + x2_traject_ij[N+1]) 
    return (point_set_x1_a_i, point_set_x1_b_i, 
            point_set_x1_a_j, point_set_x1_b_j, 
            point_set_x2_a_ij, point_set_x2_b_ij)
end;




#--------------- The followings functions are for visualization of frequency trajectory ----------#
function get_trajactory_site_i(i, q, MPLseq_raw, specific_times, muMat, pseudo_count)
    index_MPLseq = vcat(collect(1:3), km.(i,1:q,q) .+ 3)
    (temp_1, temp_2, temp_3, temp_4, x1_traject_i, x2_traject_i) = get_controle_point_for_visualization(
                q, 1, MPLseq_raw[:, index_MPLseq], specific_times, muMat, pseudo_count);
    (point_set_x1_a_B_i, point_set_x1_b_B_i, point_set_x2_a_B_i, point_set_x2_b_B_i) = 
    get_Bezier_vec_Mat_for_all_sites_time_categorical(q, 1, x1_traject_i, x2_traject_i);
    return (point_set_x1_a_B_i, point_set_x1_b_B_i, point_set_x2_a_B_i, point_set_x2_b_B_i)
end;


function get_controle_point_for_visualization(q, L, data, time_list_in, muMat, pseud_count)
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

    (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b) = 
    get_Bezier_vec_Mat_for_all_sites_time_categorical(q, L, x1_traject, x2_traject)

    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b, x1_traject, x2_traject)
end;

function check_need_to_inclue_middpoint(x1_traject, x2_traject, i, q, diff_freq_thresh=0.6)    
    #if frequency change more than 0.6 % then make a frag
    idx_t_max = size(x1_traject, 1)
    x1_traject_i = []
    x2_traject_ii = []
    
    push!(x1_traject_i, copy(x1_traject[1][km.(i,1:q,q)]))
    push!(x2_traject_ii, copy(x2_traject[1][km.(i,1:q,q), km.(i,1:q,q)]))
    for i_t in 2:idx_t_max
        flag_inserted = false
        push!(x1_traject_i, copy(x1_traject[i_t][km.(i,1:q,q)]))
        push!(x2_traject_ii, copy(x2_traject[i_t][km.(i,1:q,q), km.(i,1:q,q)]))
        for a in 1:q
            diff_freq = abs(x1_traject[i_t][km(i,a,q)] - x1_traject[i_t-1][km(i,a,q)])
            if(diff_freq>=diff_freq_thresh)
                flag_inserted = true
                break;
            end
        end
    end
    return (flag_inserted, x1_traject_i, x2_traject_ii)
end;


# suppose the i is polymorophic index.
function check_need_to_inclue_middpoint_element_wise(x1_traject, x2_traject, i, diff_freq_thresh=0.6)
    flag_inserted = false
    #if frequency change more than 0.6 % then make a frag
    idx_t_max = size(x1_traject, 1)
    x1_traject_i = []
    x2_traject_ii = []
    
    push!(x1_traject_i, copy(x1_traject[1][i]))
    push!(x2_traject_ii, copy(x2_traject[1][i, i]))
    for i_t in 2:idx_t_max
        push!(x1_traject_i, copy(x1_traject[i_t][i]))
        push!(x2_traject_ii, copy(x2_traject[i_t][i,i]))
        for a in 1:q
            diff_freq = abs(x1_traject[i_t][i] - x1_traject[i_t-1][i])
            if(diff_freq>=diff_freq_thresh)
                flag_inserted = true
                break;
            end
        end
    end
    return (flag_inserted, x1_traject_i, x2_traject_ii)
end;



function get_midpoint_freq(q, x1_traject_i, x2_traject_ii, unique_collecting_date, diff_freq_thresh = 0.6)
    #diff_freq_thresh = 0.6 #if frequency change more than 0.6 % then make a frag
    idx_t_max = size(x1_traject_i, 1)
    midpoint_vec = []
    x1_traject_i_with_midpts = []
    x2_traject_ii_with_midpts = []
    unique_collecting_date_with_midpts = []
    
    push!(x1_traject_i_with_midpts, copy(x1_traject_i[1]))
    push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[1]))
    push!(unique_collecting_date_with_midpts, unique_collecting_date[1])
    for i_t in 2:idx_t_max
        flag_inserted = false
        for a in 1:q
            diff_freq = abs(x1_traject_i[i_t][a] - x1_traject_i[i_t-1][a])
            if(diff_freq>=diff_freq_thresh && !flag_inserted)
                @show diff_freq
                midpoint_vec = 0.5*(x1_traject_i[i_t] + x1_traject_i[i_t-1])
                midpoint_mat = 0.5*(x2_traject_ii[i_t] + x2_traject_ii[i_t-1])
                midpoint_time = 0.5*(unique_collecting_date[i_t] + unique_collecting_date[i_t-1])
                push!(x1_traject_i_with_midpts, copy(midpoint_vec))
                push!(x2_traject_ii_with_midpts, copy(midpoint_mat))
                push!(unique_collecting_date_with_midpts, midpoint_time)
                flag_inserted = true
            end
        end
        push!(x1_traject_i_with_midpts, copy(x1_traject_i[i_t]))
        push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[i_t]))
        push!(unique_collecting_date_with_midpts, unique_collecting_date[i_t])
    end
    return (x1_traject_i_with_midpts, x2_traject_ii_with_midpts, unique_collecting_date_with_midpts)
end;


function get_midpoint_freq_element_wise(x1_traject_i, x2_traject_ii, unique_collecting_date, diff_freq_thresh = 0.6)
    flag_inserted = false
    idx_t_max = size(x1_traject_i, 1)
    midpoint_vec = []
    x1_traject_i_with_midpts = []
    x2_traject_ii_with_midpts = []
    unique_collecting_date_with_midpts = []
    
    push!(x1_traject_i_with_midpts, copy(x1_traject_i[1]))
    push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[1]))
    push!(unique_collecting_date_with_midpts, unique_collecting_date[1])
    
    for i_t in 2:idx_t_max
        diff_freq = abs(x1_traject_i[i_t] - x1_traject_i[i_t-1])
        if(diff_freq>=diff_freq_thresh)
            # these are not vector but scalar variable.
            midpoint_i = 0.5*(x1_traject_i[i_t] + x1_traject_i[i_t-1])
            midpoint_ii = 0.5*(x2_traject_ii[i_t] + x2_traject_ii[i_t-1])
            midpoint_time = 0.5*(unique_collecting_date[i_t] + unique_collecting_date[i_t-1])

            push!(x1_traject_i_with_midpts, copy(midpoint_i))
            push!(x2_traject_ii_with_midpts, copy(midpoint_ii))
            push!(unique_collecting_date_with_midpts, copy(midpoint_time))
            flag_inserted = true
        end
        
        push!(x1_traject_i_with_midpts, copy(x1_traject_i[i_t]))
        push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[i_t]))
        push!(unique_collecting_date_with_midpts, copy(unique_collecting_date[i_t]))
    end
    return (x1_traject_i_with_midpts, x2_traject_ii_with_midpts, unique_collecting_date_with_midpts)
end;


function get_midpoint_freq_include_neighbor(q, x1_traject_i, x2_traject_ii, unique_collecting_date, diff_freq_thresh = 0.6)
    #diff_freq_thresh = 0.6 #if frequency change more than 0.6 % then make a frag
    idx_t_max = size(x1_traject_i, 1)
    midpoint_vec = []
    x1_traject_i_with_midpts = []
    x2_traject_ii_with_midpts = []
    unique_collecting_date_with_midpts = []
    
    push!(x1_traject_i_with_midpts, copy(x1_traject_i[1]))
    push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[1]))
    push!(unique_collecting_date_with_midpts, unique_collecting_date[1])
    
    mixing_rate = 0.99 # this is rato for the actually observed frequency
    
    for i_t in 2:idx_t_max
        flag_inserted = false
        for a in 1:q
            diff_freq = abs(x1_traject_i[i_t][a] - x1_traject_i[i_t-1][a])
            if(diff_freq>=diff_freq_thresh && !flag_inserted)
                @show diff_freq
                midpoint_vec = 0.5*(x1_traject_i[i_t] + x1_traject_i[i_t-1])
                midpoint_mat = 0.5*(x2_traject_ii[i_t] + x2_traject_ii[i_t-1])
                midpoint_time = 0.5*(unique_collecting_date[i_t] + unique_collecting_date[i_t-1])
                push!(x1_traject_i_with_midpts, copy(midpoint_vec))
                push!(x2_traject_ii_with_midpts, copy(midpoint_mat))
                push!(unique_collecting_date_with_midpts, midpoint_time)
                flag_inserted = true
                
                #insert midpoint before this interval if it defined more than 2-3th points
                if(i_t>=3)
                    midpoint_vec = mixing_rate*x1_traject_i[i_t-1] + (1-mixing_rate)*x1_traject_i[i_t-2]
                    midpoint_mat = mixing_rate*x2_traject_ii[i_t-1] + (1-mixing_rate)*x2_traject_ii[i_t-2]
                    midpoint_time = mixing_rate*unique_collecting_date[i_t-1] + (1-mixing_rate)*unique_collecting_date[i_t-2]
                    push!(x1_traject_i_with_midpts, copy(midpoint_vec))
                    push!(x2_traject_ii_with_midpts, copy(midpoint_mat))
                    push!(unique_collecting_date_with_midpts, midpoint_time)
                end

                #insert midpoint after this interval if neighboring interval i_t+1 <= idx_t_max
                if(i_t+1 <= idx_t_max)
                    midpoint_vec = (1-mixing_rate)*x1_traject_i[i_t+1] + mixing_rate*x1_traject_i[i_t]
                    midpoint_mat = (1-mixing_rate)*x2_traject_ii[i_t+1] + mixing_rate*x2_traject_ii[i_t]
                    midpoint_time = (1-mixing_rate)*unique_collecting_date[i_t+1] + mixing_rate*unique_collecting_date[i_t]
                    push!(x1_traject_i_with_midpts, copy(midpoint_vec))
                    push!(x2_traject_ii_with_midpts, copy(midpoint_mat))
                    push!(unique_collecting_date_with_midpts, midpoint_time)
                end                
            end
        end
        push!(x1_traject_i_with_midpts, copy(x1_traject_i[i_t]))
        push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[i_t]))
        push!(unique_collecting_date_with_midpts, unique_collecting_date[i_t])
        
    end
    
    index_sorting = sortperm(unique_collecting_date_with_midpts)
    x1_traject_i_with_midpts = copy(x1_traject_i_with_midpts[index_sorting])
    x2_traject_ii_with_midpts = copy(x2_traject_ii_with_midpts[index_sorting])
    unique_collecting_date_with_midpts = copy(unique_collecting_date_with_midpts[index_sorting])
    
    return (x1_traject_i_with_midpts, x2_traject_ii_with_midpts, unique_collecting_date_with_midpts)
end;


function get_midpoint_freq_element_wise_include_neighbor(x1_traject_i, x2_traject_ii, unique_collecting_date, diff_freq_thresh = 0.6)
    flag_inserted = false
    idx_t_max = size(x1_traject_i, 1)
    midpoint_vec = []
    x1_traject_i_with_midpts = []
    x2_traject_ii_with_midpts = []
    unique_collecting_date_with_midpts = []
    
    push!(x1_traject_i_with_midpts, copy(x1_traject_i[1]))
    push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[1]))
    push!(unique_collecting_date_with_midpts, unique_collecting_date[1])
    
    mixing_rate = 0.99 # this is rato for the actually observed frequency
    
    for i_t in 2:idx_t_max
        diff_freq = abs(x1_traject_i[i_t] - x1_traject_i[i_t-1])
        if(diff_freq>=diff_freq_thresh)
            # these are not vector but scalar variable.
            midpoint_i = 0.5*(x1_traject_i[i_t] + x1_traject_i[i_t-1])
            midpoint_ii = 0.5*(x2_traject_ii[i_t] + x2_traject_ii[i_t-1])
            midpoint_time = 0.5*(unique_collecting_date[i_t] + unique_collecting_date[i_t-1])
            push!(x1_traject_i_with_midpts, copy(midpoint_i))
            push!(x2_traject_ii_with_midpts, copy(midpoint_ii))
            push!(unique_collecting_date_with_midpts, copy(midpoint_time))
            flag_inserted = true
            
            #insert midpoint before this interval if it defined more than 2-3th points
            if(i_t>=3)
                midpoint_i = mixing_rate*x1_traject_i[i_t-1] + (1-mixing_rate)*x1_traject_i[i_t-2]
                midpoint_ii = mixing_rate*x2_traject_ii[i_t-1] + (1-mixing_rate)*x2_traject_ii[i_t-2]
                midpoint_time = mixing_rate*unique_collecting_date[i_t-1] + (1-mixing_rate)*unique_collecting_date[i_t-2]
                push!(x1_traject_i_with_midpts, copy(midpoint_i))
                push!(x2_traject_ii_with_midpts, copy(midpoint_ii))
                push!(unique_collecting_date_with_midpts, copy(midpoint_time))                
            end
            
            #insert midpoint after this interval if neighboring interval i_t+1 <= idx_t_max
            if(i_t+1 <= idx_t_max)
                midpoint_i = (1-mixing_rate)*x1_traject_i[i_t+1] + mixing_rate*x1_traject_i[i_t]
                midpoint_ii = (1-mixing_rate)*x2_traject_ii[i_t+1] + mixing_rate*x2_traject_ii[i_t]
                midpoint_time = (1-mixing_rate)*unique_collecting_date[i_t+1] + mixing_rate*unique_collecting_date[i_t]
                push!(x1_traject_i_with_midpts, copy(midpoint_i))
                push!(x2_traject_ii_with_midpts, copy(midpoint_ii))
                push!(unique_collecting_date_with_midpts, copy(midpoint_time))                
            end
        end        
        push!(x1_traject_i_with_midpts, copy(x1_traject_i[i_t]))
        push!(x2_traject_ii_with_midpts, copy(x2_traject_ii[i_t]))
        push!(unique_collecting_date_with_midpts, copy(unique_collecting_date[i_t]))
    end
    
    index_sorting = sortperm(unique_collecting_date_with_midpts)
    x1_traject_i_with_midpts = copy(x1_traject_i_with_midpts[index_sorting])
    x2_traject_ii_with_midpts = copy(x2_traject_ii_with_midpts[index_sorting])
    unique_collecting_date_with_midpts = copy(unique_collecting_date_with_midpts[index_sorting])
    
    return (x1_traject_i_with_midpts, x2_traject_ii_with_midpts, unique_collecting_date_with_midpts)
end;
