function priority = compute_priority(W_parent, W_child)

n = length(W_parent);
[sorted_parent, idx_parent] = sort(W_parent, 'descend');
[sorted_child1, idx_child1] = sort(W_child(:, 1), 'descend');
[sorted_child2, idx_child2] = sort(W_child(:, 2), 'descend');

n_part = length(find(W_parent ~= 0));
if n_part <= 1
	priority = -3;
else
	weight = log(n:-1:1)';
	first_zero = find(sorted_parent == 0, 1); 
	if length(first_zero) > 0 
		weight(first_zero:end) = 1;
	end
	weight_part = zeros(n, 1);
	weight_part(1:n_part) = log(n_part:-1:1)';
	[sorted, idx1] = sort(idx_child1);
	[sorted, idx2] = sort(idx_child2);
	max_pos = max(idx1, idx2);
	discount = log(n-max_pos(idx_parent)+1);
	discount(discount == 0) = log(2);
	weight = weight ./ discount;
	weight_part = weight_part ./ discount;
	priority = NDCG_part(idx_parent, idx_child1, weight, weight_part) * NDCG_part(idx_parent, idx_child2, weight, weight_part);
end