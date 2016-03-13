function [cluster_subset, W_buffer_one, H_buffer_one, priority_one] = actual_split(X, subset, W_parent, params)

[m, n] = size(X);
if length(subset) <= 3
	cluster_subset = ones(1, length(subset));
	W_buffer_one = zeros(m, 2);
	H_buffer_one = zeros(2, length(subset));
	priority_one = -1;
else
	term_subset = find(sum(X(:, subset), 2) ~= 0);
	X_subset = X(term_subset, subset);
	W = rand(length(term_subset), 2);
	H = rand(2, length(subset));
	[W, H] = nmfsh_comb_rank2(X_subset, W, H, params);
	[max_val, cluster_subset] = max(H);
	W_buffer_one = zeros(m, 2);
	W_buffer_one(term_subset, :) = W;
	H_buffer_one = H;
	if length(unique(cluster_subset)) > 1
		priority_one = compute_priority(W_parent, W_buffer_one);
	else
		priority_one = -1;
	end
end