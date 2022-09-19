# This function returns samples at the time t
function get_sample_at_t(data, t_get)
    n_t = []
    sample_t = []
    for n in 1:size(data,1)
        if(Int(data[n, 1]) == t_get)
            if(length(sample_t)>0)
                sample_t = hcat(sample_t, data[n,4:end]) 
                n_t = vcat(n_t, data[n,2])
            end
            if(length(sample_t)==0)
                sample_t = copy(data[n,4:end])
                n_t = data[n,2]
            end
        end

        if(data[n, 1] > t_get)
            break
        end
    end
    return (n_t, sample_t)
end;



function get_TPR_curve_bene(s_ben, s_neu, s_dele)
    s_tot = [s_ben; s_neu; s_dele]
    s_label_beneficial = [ones(size(s_ben)); zeros(size(s_neu)); zeros(size(s_dele))];
    s_label_deleterious = [zeros(size(s_ben)); zeros(size(s_neu)); ones(size(s_dele))];
    index_TP = sortperm(s_tot, rev=true);
    index_FP = sortperm(s_tot);

    #Note FN = 1-TP, TN = 1-FP 
    tp_bene = s_label_beneficial[index_TP] # Beneficial TP
    #fp_bene = s_label_beneficial[index_FP] # Beneficial FP;
    tp_dele = s_label_deleterious[index_FP] # Deleterious TP
    #fp_dele = s_label_deleterious[index_TP]; # Deleterious FP;

    size_tot_length = size(tp_bene,1)
    TPR_bene = zeros(size_tot_length)
    TPR_dele = zeros(size_tot_length)

    count_TP_bene = 0 
    count_TP_dele = 0 

    True_bene_tot = size(s_ben,1)
    True_dele_tot = size(s_dele,1)

    for n in 1:size(TPR_bene,1)
        count_TP_bene += tp_bene[n]
        TPR_bene[n] = count_TP_bene*(1.0/n)

        count_TP_dele += tp_dele[n]
        TPR_dele[n] = count_TP_dele*(1.0/n)
    end
    return (TPR_bene, TPR_dele)
end

linreg(x,y) = (hcat(fill!(similar(x),1),x) \y)[2]; # This gives different result from the minimization of square error.
linreg_expricit(x,y) = (x'*x)\x'*y;
linreg_tot(x,y) = (hcat(fill!(similar(x),1),x) \y) ;

# Functions for obtaining PPV curves for WF-MPL
function read_coupling(fname_in)
    J_raw = readdlm(fname_in)
    L = Int(J_raw[end,1])
    J_in = zeros(L, L)
    for n in 1:size(J_raw, 1)
        x = J_raw[n, :]
        i,j,v = Int(x[1]), Int(x[2]), x[3]
        J_in[i,j] = v
        J_in[j,i] = v
    end
    return J_in
end
function get_diag_offdiag(J)
    J_diag = J[diagind(J)]
    J_offdiag = J - diagm(0=>J_diag)
    return(J_diag, J_offdiag)
end;

function get_corslope_label(vec_J1, vec_J2)
    #slope = linreg(vec_J1, vec_J2)
    slope = linreg_expricit(vec_J1, vec_J2)
    cor = Statistics.cor(vec_J1, vec_J2)
    return text_string = @sprintf("Slope=%.3f, Corr.=%.3f\n", slope, cor)
end;

function get_corslope_label_v2(vec_J1, vec_J2)
    slope = linreg(vec_J1, vec_J2)
    #cor = Statistics.cor(vec_J1, vec_J2)
    return text_string = @sprintf("Slope=%.3f\n", slope)
end;


function get_corslope(vec_J1, vec_J2)
    slope = linreg(vec_J1, vec_J2)
    cor = Statistics.cor(vec_J1, vec_J2)
    return (slope, cor)
end;
