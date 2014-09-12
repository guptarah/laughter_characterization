function [all_feats] = get_code_vectors_around_laugh_v3(win_size,num_feats)

% In this version we get last n counselor and client codes for these 5 classes
% 1: no laughter
% 2: stand alone client laughter
% 3: stand alone counselor laughter
% 4: client laughter followed by couns laughter
% 5: couns laughter followed by client laughter  
% all_feats : [file_id class features]
files01 = dir('../../data/to_use_set_code01/*.turn');

%win_size = 5;
couns_all_data = [];
client_all_data = [];

client_to_couns_dist = [];
couns_to_client_dist = [];
shared_laugh_feats = [];
stand_alone_laugh_feats = [];
all_feats = [];
for file = files01'
	
	cur_file=strcat('../../data/to_use_set_code01/',file.name);
	disp(file.name);
	cur_data = load (cur_file);
	file_id = str2num(strtok(file.name,'_'));

	ext_cur_data = [zeros(win_size,3); cur_data; zeros(1,3)];
	
	for i = 1:size(cur_data,1)
		cur_feats = ext_cur_data(i:i+win_size,1);
		cur_feats1 = cur_feats(find(and((cur_feats > 10),(cur_feats<20))));
		cur_feats2 = cur_feats(find(cur_feats>20));

		if length(cur_feats1) >= num_feats
			cur_feats1 = cur_feats1(end-num_feats+1:end);
		else
			cur_feats1 = [zeros(num_feats-length(cur_feats1),1) ; cur_feats1];
		end

		if length(cur_feats2) >= num_feats
                        cur_feats2 = cur_feats2(end-num_feats+1:end);
                else
                        cur_feats2 = [zeros(num_feats-length(cur_feats2),1) ; cur_feats2];
                end

		cur_feats = [cur_feats1;cur_feats2];

		cur_has_laugh = ext_cur_data(i+win_size,3);
		cur_speaker = floor(ext_cur_data(i+win_size,1)/20); % 0 meaning couns 1 meaninig client		
		next_speaker = floor(ext_cur_data(i+win_size+1,1)/20);
		next_has_laugh = ext_cur_data(i+win_size+1,3);
	
		% assigning classes
		if cur_has_laugh == 0
			class = 1;
		elseif cur_speaker == 0 % counselor current 
			if (next_speaker == 1) && (next_has_laugh == 1)
				class = 5;
			else
				class = 3;	
			end
		elseif cur_speaker == 1 % client current
			if (next_speaker == 0) && (next_has_laugh == 1)
                                class = 4;
                        else
                                class = 2;
                        end
		end
	
		all_feats = [all_feats; file_id class cur_feats']; 	
	end
	
end

%bigram_all_feats = compute_bigrams(all_feats);
csvwrite('all_feats.v3',all_feats);
