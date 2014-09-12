function [all_feats] = get_code_vectors_around_laugh_v4(win_size)

% Features are the counts of each code
% In this version we get for these 5 classes
% 1: no laughter
% 2: stand alone client laughter
% 3: stand alone counselor laughter
% 4: client laughter followed by couns laughter
% 5: couns laughter followed by client laughter  
% all_feats : [file_id class features]
files01 = dir('../../data/to_use_set_code01/*.turn');
delete all_feats.v4;

%win_size = 5;
couns_all_data = [];
client_all_data = [];

client_to_couns_dist = [];
couns_to_client_dist = [];
shared_laugh_feats = [];
stand_alone_laugh_feats = [];
all_feats = [];

code_list = [11:15 21:23];
bigram_list = repmat(100*code_list',1,8)+repmat(code_list,8,1);
bigram_list = reshape(bigram_list,1,64);

for file = files01'
	
	cur_file=strcat('../../data/to_use_set_code01/',file.name);
	disp(file.name);
	cur_data = load (cur_file);
	file_id = str2num(strtok(file.name,'_'));

	ext_cur_data = [zeros(win_size,3); cur_data; zeros(1,3)];
	
	for i = 1:size(cur_data,1)
		cur_feats = ext_cur_data(i:i+win_size,1);
		bigrams = 100*cur_feats(1:end-1) + cur_feats(2:end);

		mono_feat = zeros(1,8);
		for code_id = 1:8 
			%temp_feat(code_id) = code_list(code_id)+sum(cur_feats == code_list(code_id))/(win_size+1);
			mono_feat(code_id) = code_list(code_id)+sum(cur_feats == code_list(code_id))/10;	
		end
		
		
		bi_feat = zeros(1,8^2);
		for code_id = 1:8^2
			bi_feat(code_id) = bigram_list(code_id) + sum(bigrams == bigram_list(code_id))/10;
		end	

		cur_feats = [mono_feat bi_feat];

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
	
		all_feats = [all_feats; file_id class cur_feats]; 
		%dlmwrite('all_feats.v4',all_feats,'-append','delimiter',',','precision','%.1f');
	end
	
	dlmwrite('all_feats.v4',all_feats,'-append','delimiter',',','precision','%.1f');
	all_feats = [];	
end

%dlmwrite('all_feats.v4',all_feats,'-append','delimiter',',','precision','%.1f');

%bigram_all_feats = all_feats;
%bigram_all_feats = compute_bigrams(all_feats);
%csvwrite('all_feats.v4',bigram_all_feats);
%dlmwrite('all_feats.v4',bigram_all_feats,'delimiter',',','precision','%.1f');
