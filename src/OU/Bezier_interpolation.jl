## This matrix is for the x_Bez and d/dt x_Bez product
get_integrated_Bezier_single() = (3.0/60) * [
    -10 -2 -1 -1;
    2 0 -1/3 -1;
    1 1/3 0 -2;
    1 1 2 10] ;

function get_Bezier_vec_Mat_for_all_sites_time(L, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    point_set_x1 = zeros(L, n_traject)
    point_set_x1_a = zeros(L, N)
    point_set_x1_b = zeros(L, N)
    for i in 1:L
        point_set_x1[i,:] = [ x1_traject[n][i] for n in 1:n_traject];
        size(point_set_x1), size(point_set_x1[i,:]), N
        (a_vec_i, b_vec_i) = get_a_b_Bezier(N, point_set_x1[i,:])
        point_set_x1_a[i,:] = copy(a_vec_i)
        point_set_x1_b[i,:] = copy(b_vec_i)
    end
    #--------------------------#
    #-- Interpolation points for corr. for each time --#
    point_set_x2 = zeros(L,L, n_traject)
    point_set_x2_a = zeros(L,L, N)
    point_set_x2_b = zeros(L,L, N)

    for i in 1:L
        for j in i:L
            point_set_x2[i,j,:] = [ x2_traject[n][i,j] for n in 1:n_traject];
            (a_vec_i, b_vec_i) = get_a_b_Bezier(N, point_set_x2[i,j,:])
            point_set_x2_a[i,j,:] = copy(a_vec_i)
            point_set_x2_a[j,i,:] = copy(a_vec_i)
            point_set_x2_b[j,i,:] = copy(b_vec_i)
            point_set_x2_b[i,j,:] = copy(b_vec_i)

        end
    end
    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;

function get_Bezier_vec_Mat_for_all_sites_time_2ndOrder(L, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    point_set_x1 = zeros(L, n_traject)
    point_set_x1_a = zeros(L, N)
    for i in 1:L
        point_set_x1[i,:] = [ x1_traject[n][i] for n in 1:n_traject];
        size(point_set_x1), size(point_set_x1[i,:]), N
        (a_vec_i) = get_a_Bezier_2ndOrder(point_set_x1[i,:])
        point_set_x1_a[i,:] = copy(a_vec_i)
    end
    #--------------------------#
    #-- Interpolation points for corr. for each time --#
    point_set_x2 = zeros(L,L, n_traject)
    point_set_x2_a = zeros(L,L, N)

    for i in 1:L
        for j in i:L
            point_set_x2[i,j,:] = [ x2_traject[n][i,j] for n in 1:n_traject];
            (a_vec_i) = get_a_Bezier_2ndOrder(point_set_x2[i,j,:])
            point_set_x2_a[i,j,:] = copy(a_vec_i)
            point_set_x2_a[j,i,:] = copy(a_vec_i)

        end
    end
    return (point_set_x1_a, point_set_x2_a)
end;



"""
Retuns the interpolation: n_out is the number of points that can interpolate between the observation points [p_n, p_{n+1}]. 
"""
function get_bezier_interpolation_periodic(n_out, x_set, p_set_left, p_set_right, a_vec, b_vec)
    scale_n_out = 1.0 / (n_out-1)
    x_axis = collect(x_set[1]:(x_set[2]-1))
    y_axis = Bez.( scale_n_out*collect(0:(n_out-1)), p_set_left[1], p_set_right[1], a_vec[1], b_vec[1]);
    for n in 1:(N-1)
        x_axis = vcat(x_axis, collect(x_set[n+1]:(x_set[n+2]-1)) )
        y_axis = vcat(y_axis, Bez.( scale_n_out*collect(0:(n_out-1)), p_set_left[n+1], p_set_right[n+1], a_vec[n+1], b_vec[n+1]))
    end
    return (x_axis, y_axis)
end
