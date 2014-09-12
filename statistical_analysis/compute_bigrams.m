function all_feats = compute_bigrams(all_feats)

%all_feats = load(all_feat_file);

features = all_feats(:,3:end);

first_part = features(:,1:end-1);
second_part = features(:,2:end);

bigrams = 100*first_part + second_part;

all_feats = [all_feats bigrams];
%csvwrite('all_feats.v2_bigrams',all_feats);
