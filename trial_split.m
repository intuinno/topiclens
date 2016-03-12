function [subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, X, subset, W_parent, params)

[m, n] = size(X);

trial = 0;
subset_backup = subset;
while trial < trial_allowance
	[cluster_subset, W_buffer_one, H_buffer_one, priority_one] = actual_split(X, subset, W_parent, params);
	if priority_one < 0
		break;
	end
	unique_cluster_subset = unique(cluster_subset);
	if length(unique_cluster_subset) ~= 2
		error('Invalid number of unique sub-clusters!');
	end
	length_cluster1 = length(find(cluster_subset == unique_cluster_subset(1)));
	length_cluster2 = length(find(cluster_subset == unique_cluster_subset(2)));
	if min(length_cluster1, length_cluster2) < unbalanced * length(cluster_subset)
		[min_val, idx_small] = min([length_cluster1, length_cluster2]);
		subset_small = find(cluster_subset == unique_cluster_subset(idx_small));
		subset_small = subset(subset_small);
		[cluster_subset_small, W_buffer_one_small, H_buffer_one_small, priority_one_small] = actual_split(X, subset_small, W_buffer_one(:, idx_small), params);
		if priority_one_small < min_priority
			trial = trial + 1;
			if trial < trial_allowance
% 				disp(['Drop ', num2str(length(subset_small)), ' documents ...']);
				subset = setdiff(subset, subset_small);
			end
		else
			break;
		end
	else
		break;
	end
end

if trial == trial_allowance
% 	disp(['Recycle ', num2str(length(subset_backup) - length(subset)), ' documents ...']);
	subset = subset_backup;
	W_buffer_one = zeros(m, 2);
	H_buffer_one = zeros(2, length(subset));
	priority_one = -2;
end