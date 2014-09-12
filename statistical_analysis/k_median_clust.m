function [clust_medians,asng] =  k_median_clust(data,K)

% data: data to cluster
% K: number of clusters

N = size(data,1);
D = size(data,2);

% Randomly assign data to clusters
asng = randi(K,N,1);
clust_medians = zeros(K,D);

counter = 0;
sum_dist = [];
flag = 0;
while flag == 0

	counter = counter+1;
	% compute the median of each cluster
	for k = 1:K
		cur_data = data(asng == k,:);
		clust_medians(k,:) = median(cur_data);
	end

	% compute distances from new cluster points
	distances = 999*ones(N,K);
	for k = 1:K
		cur_med = clust_medians(k,:);
		% 0-1 cost function
		% distances(:,k) = sum(((data - repmat(cur_med,N,1))~=0),2);
		% L1 distance
		distances(:,k) = sum(abs(data - repmat(cur_med,N,1)),2);
	end

	% compute new asng
	[min_distances,asng] = min(distances');	
	asng = asng';
	sum_dist = [sum_dist sum(min_distances)];
	disp(sum(min_distances));
	
	if counter > 100
		flag=1;
	end

end

plot(sum_dist,'*');
