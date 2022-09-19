
function get_norm_measures(C_tot)
    norms_vec = zeros(6)
    norms_vec[1] = norm(C_tot[diagind(C_tot)][1:10]) #norm_diag_B
    norms_vec[2] = norm(C_tot[diagind(C_tot)][11:40]) #norm_diag_N
    norms_vec[3] = norm(C_tot[diagind(C_tot)][41:end]) #norm_diag_D
    norms_vec[4] = norm(C_tot - diagm( 0=>C_tot[diagind(C_tot)] ) ) #norm_offdiag
    norms_vec[5] = LinearAlgebra.tr(C_tot) #norm_trace
    norms_vec[6] = det(C_tot+0.1*I) #norm_logdet
    return norms_vec
end;


