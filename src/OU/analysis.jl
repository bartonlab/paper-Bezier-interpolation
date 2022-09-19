
function get_patterns_J(L,P)
    patterns_embbeded = 2*rand(0:1, L, P) .- 1 + 0.3*rand(L, P);
    J = - patterns_embbeded * patterns_embbeded' / (L*sqrt(P)) ;
    return (patterns_embbeded, J)
end;

get_energy(x,J) = -0.5 * x' * J * x;

function get_diag_offdiag(J)
    J_diag = J[diagind(J)]
    J_offdiag = J - diagm(0=>J_diag)
    return(J_diag, J_offdiag)
end;

function get_corslope_label(vec_J1, vec_J2)
    slope = linreg(vec_J1, vec_J2)
    cor = Statistics.cor(vec_J1, vec_J2)
    return text_string = @sprintf("slope=%.3f, corr.=%.3f\n", slope, cor)
end;

function get_corslope(vec_J1, vec_J2)
    slope = linreg(vec_J1, vec_J2)
    cor = Statistics.cor(vec_J1, vec_J2)
    return (slope, cor)
end;

function get_x_t_ensembles(L, T, dt, D, n_ensemble, J)
    sqrt_Ddt = sqrt(dt*D)
    t_ensembles = []
    x_ensembles = zeros(n_ensemble, T+1, L)

    for n in 1:n_ensemble
        x0 = 10*rand(-0.5:0.5, L);
        x = copy(x0)
        for n_t in 0:T
            x += dt * J*copy(x) + sqrt_Ddt* randn(L)
            x_ensembles[n, n_t+1, :] = copy(real(x))
            if(n==1)
                t_ensembles = push!(t_ensembles, n_t*dt)
            end
        end
    end
    return (x_ensembles, t_ensembles)
end

function get_effective_energy_OU(n_ensemble, n_t_observe, T, L, J, sqrt_Ddt)

    x_ensembles_temp = zeros(n_ensemble, T, L)
    t_ensembles_temp = []

    x0 = zeros(L)
    for i in 1:(Int(floor(0.3*L)))
        x0[i] = 5
    end
    for i in (Int(floor(0.7*L))):L
        x0[i] = -5
    end

    n_t_count_max = 1
    for n in 1:n_ensemble
        x = copy(x0)
        n_t_count = 1
        for n_t in 0:T
            #x += dt * ( J*copy(x) + sqrt_Ddt*randn(L) )  
            x += dt * J*copy(x) +sqrt_Ddt* randn(L) 
            if(n_t%n_t_observe==0)
                x_ensembles_temp[n, n_t_count, :] = copy(x)
                if(n==1)
                    t_ensembles_temp = push!(t_ensembles_temp, n_t * dt)
                end
                n_t_count +=1
            end
        end
        if(n==1)
            n_t_count_max = n_t_count-1
        end
    end
    return erg_evolution_avg = [mean([get_energy(x_ensembles_temp[n, n_t_count,:], J) for n in 1:n_ensemble_temp]) for n_t_count in 1:n_t_count_max];
end;

function estimate_Linear(L, x_ensembles, t_ensembles)
    my_numerator_Lin = zeros(L,L)
    my_denominator_Lin = zeros(L,L);
    n_ensemble = size(x_ensembles,1)
    K = size(x_ensembles,2);
    for n in 1:n_ensemble
        #for k in 1:(K_cut-1)
        for k in 1:(K-1)
            dt_int = t_ensembles[k+1] - t_ensembles[k]

            x_k = copy(x_ensembles[n,k,:])
            x_kp1 = copy(x_ensembles[n,k+1,:])

            ##### integrated outer product of xΔx^T is (x(1)+x(0))*(x(1)-x(0))^T
            my_numerator_Lin += copy(0.5 * (x_kp1+x_k)*(x_kp1-x_k)' ) 
            ##### integrated covariance #####
            my_denominator_Lin += copy(dt_int * x_ensembles[n, k:(k+1),:]' * Q1 * x_ensembles[n, k:(k+1),:] )
        end
    end
    return J_est_Lin = (my_denominator_Lin + I) \  my_numerator_Lin ;
end;


function estimate_Bezier(L, x_ensembles, t_ensembles)
    my_numerator_Bez = zeros(L,L)
    my_denominator_Bez = zeros(L,L);
    n_ensemble = size(x_ensembles,1)

    K = size(x_ensembles,2);
    for n in 1:n_ensemble
        #estimating the controle point vectos a and b.
        a_vec_set = zeros(K-1, L)
        b_vec_set = zeros(K-1, L)
        for i in 1:L
            (a_vec, b_vec) = get_a_b_Bezier(K-1, x_ensembles[n,:,i])
            a_vec_set[:, i] = copy(a_vec)
            b_vec_set[:, i] = copy(b_vec)    
        end

        #for k in 1:(K_cut-1)
        for k in 1:(K-1)            
            dt_int = t_ensembles[k+1] - t_ensembles[k]

            #----- set x_Bez matrix ------#
            x_Bez = zeros(L, 4)
            x_Bez[:, 1] = copy(x_ensembles[n,k,:])
            x_Bez[:, 2] = 3*copy(a_vec_set[k, :]) # coefficient 3 should be carefully treated!
            x_Bez[:, 3] = 3*copy(b_vec_set[k, :])
            x_Bez[:, 4] = copy(x_ensembles[n,k+1,:]);

            #------- x_Bez and x_Bez product for the denominator -------#
            x_xT_Bez = x_Bez * Mat_Bez * x_Bez';
            #------- x_Bez and d/dt x_Bez product for the numerator -------#
            x_dxT_Bez = (x_Bez * Mat_Bez_2 * x_Bez')';

            ##### integrated outer product of xΔx^T is (x(1)+x(0))*(x(1)-x(0))^T
            my_numerator_Bez += copy( x_dxT_Bez ) 
            ##### integrated covariance #####
            my_denominator_Bez += copy( dt_int * x_xT_Bez )
        end
    end
    J_est_Bez = (my_denominator_Bez + I) \  my_numerator_Bez ;
end;

