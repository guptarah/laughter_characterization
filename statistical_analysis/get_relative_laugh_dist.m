function get_relative_laugh_dist()

files01 = dir('../../data/to_use_set_01/*.turn');

client_to_couns_dist = [];
couns_to_client_dist = [];
num_cl_laughs = 0;
num_cl_empty = 0;
num_couns_laughs = 0;
for file = files01'
	
	cur_file=strcat('../../data/to_use_set_01/',file.name);
	disp(file.name);
	cur_data = load (cur_file);
	
	counselor_laugh_loc = find((cur_data(:,1) == 0).*cur_data(:,3));
	client_laugh_loc = find((cur_data(:,1) == 1).* cur_data(:,3));

%	num_cl_laughs = num_cl_laughs + length(client_laugh_loc);
	num_couns_laughs = num_couns_laughs + length(counselor_laugh_loc);

	% computing closest couns laugh to client ones
	for i = 1:length(client_laugh_loc)
		num_cl_laughs = num_cl_laughs + 1; 
		cur_laugh_loc = client_laugh_loc(i);
		compare_mat = ones(size(counselor_laugh_loc))*cur_laugh_loc;
		if length(compare_mat) == 0
			num_cl_empty = num_cl_empty + 1;
		end 
		laugh_distances = counselor_laugh_loc - compare_mat;
                [~,closest_laugh_loc] = min(abs(laugh_distances));
                closest_laugh = laugh_distances(closest_laugh_loc);
%		disp(closest_laugh);
		couns_to_client_dist = [couns_to_client_dist closest_laugh];
	end

	% computing closest client laugh to couns ones
	for i = 1:length(counselor_laugh_loc)

                cur_laugh_loc = counselor_laugh_loc(i);
                compare_mat = ones(size(client_laugh_loc))*cur_laugh_loc;
		laugh_distances = client_laugh_loc - compare_mat;
                [~,closest_laugh_loc] = min(abs(laugh_distances));
		closest_laugh = laugh_distances(closest_laugh_loc);
%                disp(closest_laugh);
                client_to_couns_dist = [client_to_couns_dist closest_laugh];
        end

end
