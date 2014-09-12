function [couns_all_data,client_all_data] = get_code_vectors_around_laugh()

# In this version we get client/couns shared vs stand alone laughters with codes in surrounding window

files01 = dir('../../data/to_use_set_code01/*.turn');

win_size = 5;
couns_all_data = [];
client_all_data = [];

client_to_couns_dist = [];
couns_to_client_dist = [];
shared_laugh_feats = [];
stand_alone_laugh_feats = [];
for file = files01'
	
	cur_file=strcat('../../data/to_use_set_code01/',file.name);
	disp(file.name);
	cur_data = load (cur_file);
	
	counselor_laugh_loc = find((floor(cur_data(:,1)/10) == 1).*cur_data(:,3));
	client_laugh_loc = find((floor(cur_data(:,1)/20) == 1).* cur_data(:,3));

	% computing client features 
	for i = 1:length(client_laugh_loc)
		cur_laugh_loc = client_laugh_loc(i);
		% now finding if laughter is shared or not
		compare_mat = ones(size(counselor_laugh_loc))*cur_laugh_loc;
		laugh_distances = counselor_laugh_loc - compare_mat;
                [min_dist,closest_laugh_loc] = min(abs(laugh_distances));
                
		if min_dist == 1 % meaning the laughter is shared
			file_id = str2num(strtok(file.name,'_'));
			laugh_person = 1;
			shared = 1;
			
			% taking care if laughter occurs at very end of begining
			if ~(cur_laugh_loc-win_size < 1 || cur_laugh_loc+win_size > size(cur_data,1))
				features = cur_data(cur_laugh_loc-win_size:cur_laugh_loc+win_size,1)';
			elseif cur_laugh_loc-win_size < 1 
				features = cur_data(1:cur_laugh_loc+win_size,1)';
				len_zeros = 2*win_size+1 - length(features);
				features = [zeros(1,len_zeros) features];
			elseif cur_laugh_loc+win_size > size(cur_data,1)
				features = cur_data(cur_laugh_loc-win_size:end,1)';
				len_zeros = 2*win_size+1 -length(features);
				features = [features zeros(1,len_zeros)]; 
			end

			client_all_data = [client_all_data; file_id laugh_person shared features]; 
		else
			file_id = str2num(strtok(file.name,'_'));
			laugh_person = 1;
			shared = 0;
			
			if ~(cur_laugh_loc-win_size < 1 || cur_laugh_loc+win_size > size(cur_data,1))
                                features = cur_data(cur_laugh_loc-win_size:cur_laugh_loc+win_size,1)';
                        elseif cur_laugh_loc-win_size < 1
                                features = cur_data(1:cur_laugh_loc+win_size,1)';
                                len_zeros = 2*win_size+1 - length(features);
                                features = [zeros(1,len_zeros) features];
                        elseif cur_laugh_loc+win_size > size(cur_data,1)
                                features = cur_data(cur_laugh_loc-win_size:end,1)';
                                len_zeros = 2*win_size+1 -length(features);
                                features = [features zeros(1,len_zeros)];
                        end	
			
			client_all_data = [client_all_data; file_id laugh_person shared features];
		end
	end

	% computing  
	for i = 1:length(counselor_laugh_loc)

                cur_laugh_loc = counselor_laugh_loc(i);
                compare_mat = ones(size(client_laugh_loc))*cur_laugh_loc;
		laugh_distances = client_laugh_loc - compare_mat;
                [min_dist,closest_laugh_loc] = min(abs(laugh_distances));

		if min_dist == 1 % meaning the laughter is shared
                        file_id = str2num(strtok(file.name,'_'));
                        laugh_person = 0;
                        shared = 1;
                        
			if ~(cur_laugh_loc-win_size < 1 || cur_laugh_loc+win_size > size(cur_data,1))
                                features = cur_data(cur_laugh_loc-win_size:cur_laugh_loc+win_size,1)';
                        elseif cur_laugh_loc-win_size < 1
                                features = cur_data(1:cur_laugh_loc+win_size,1)';
                                len_zeros = 2*win_size+1 - length(features);
                                features = [zeros(1,len_zeros) features];
                        elseif cur_laugh_loc+win_size > size(cur_data,1)
                                features = cur_data(cur_laugh_loc-win_size:end,1)';
                                len_zeros = 2*win_size+1 -length(features);
                                features = [features zeros(1,len_zeros)];
                        end

			couns_all_data = [couns_all_data; file_id laugh_person shared features];
                else
                        file_id = str2num(strtok(file.name,'_'));
                        laugh_person = 0;
                        shared = 0;
                        
			if ~(cur_laugh_loc-win_size < 1 || cur_laugh_loc+win_size > size(cur_data,1))
                                features = cur_data(cur_laugh_loc-win_size:cur_laugh_loc+win_size,1)';
                        elseif cur_laugh_loc-win_size < 1 
                                features = cur_data(1:cur_laugh_loc+win_size,1)';
                                len_zeros = 2*win_size+1 - length(features);
                                features = [zeros(1,len_zeros) features];
                        elseif cur_laugh_loc+win_size > size(cur_data,1)
                                features = cur_data(cur_laugh_loc-win_size:end,1)';
                                len_zeros = 2*win_size+1 -length(features);
                                features = [features zeros(1,len_zeros)];
                        end

			couns_all_data = [couns_all_data; file_id laugh_person shared features];
                end
        end

end
