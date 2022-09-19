# MPL Tools
get_time_list(data) = Int.(unique(data[:,1]));

# ------ this function enhacne the genetic drift ---- 
#N_subset is a number of the populatin size which is smaller than original population size N.
function resumpling_data(data, N, N_subset)
	data_resumpled = zeros(size(data));
	n_shift = 0
	for t_sampling in unique(data[:, 1])
	    #----- extract a specific time dependent subsamples
	    data_temp = data[ collect(1:size(data,1))[data[:,1] .== t_sampling],  : ]
	    #----- extract the number and id of sequences 
	    temp_index = []
	    for n in 1:size(data_temp,1)
		temp_index = vcat(temp_index, [n for _ in 1:Int(data_temp[n,2]) ])
	    end
	    #------ Selecting the subsamples
	    index_subsamples = temp_index[shuffle(collect(1:N))[1:N_subset]];
	    #------ records the sequence id in the t-dependent subset 
	    index_subsamples_unique = unique(index_subsamples)
	    #------ records the number of the corresponding sequences
	    n_count_subsamples = [count(index_subsamples .== n) for n in index_subsamples_unique];

	    #------ slect subsamples from the original time depending subsamples.
	    data_temp2 = data_temp[index_subsamples_unique, :]
	    #------ filling the number of the sequences
	    for n in 1:size(n_count_subsamples,1)
		data_temp2[n, 2] = n_count_subsamples[n]
	    end

	    #------ Insert resampled sequences to the output matrx.
	    data_resumpled[(n_shift+1):(n_shift+size(data_temp2,1)), :] = copy(data_temp2)
	    n_shift += size(data_temp2,1)
	end

	return data_resumpled = copy(data_resumpled[1:n_shift, :]);
end;

function get_data_time(file_key1, file_key2, data_dir, id_ensemble, time_upto=300)
    fname_in = data_dir * file_key1 *"_id-"*string(id_ensemble) * file_key2
    data = readdlm(fname_in);
    read_upto = count(data[:,1] .<= time_upto)
    data = copy(data[1:read_upto, :])
    time_list = get_time_list(data);
    return (data, time_list)
end;



