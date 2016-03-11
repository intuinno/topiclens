function [tree, splits, is_leaf, clusters, timings, Ws, priorities, W, H] = hier8_neat_reduced(X, k, cl_idx, idx, params)
%%HIER8_NEAT - Hierarchical clustering based on rank-2 NMF
% [tree, splits, is_leaf, clusters, timings, Ws, priorities] = hier8_neat(X, k, params)
% Input parameters --
% X: m*n data matrix (m features x n data points)
% k: The max number of leaf nodes to be generated
% params (optional)
% params.trial_allowance (default=3): Number of trials allowed for removing
%				      outliers and splitting a node again.
%				      See parameter T in Algorithm 3 in the reference paper.
% params.unbalanced (default=0.1): A threshold to determine if one of the two clusters is an outlier set.
%				   A smaller value means more tolerance for unbalance between two clusters.
%				   See parameter beta in Algorithm 3 in the reference paper.
% params.vec_norm (default=2): Indicates which norm to use for the normalization of W or H,
%			       e.g. vec_norm=2 means Euclidean norm; vec_norm=0 means no normalization.
% params.normW (default=true): true if normalizing columns of W; false if normalizing rows of H.
% params.anls_alg (default=@anls_entry_rank2_precompute): The function handle to NNLS algorithm
%                                                         whose signature is the same as @anls_entry_rank2_precompute
% params.tol (default=1e-4): Tolerance parameter for stopping criterion in each run of NMF.
% params.maxiter (default=10000): Maximum number of iteration times in each run of NMF.
%
% Output parameters --
% From the output parameters, you can reconstruct the tree and "replay" the k-1 steps that generated it.
%
% For a binary tree with k leaf nodes, the total number of nodes (including leaf and non-leaf nodes)
% is 2*(k-1) plus the root node, because k-1 splits are performed and each split generates two new nodes.
%
% We only keep track of the 2*(k-1) non-root node in the output.
%
% tree: A 2-by-(k-1) matrix that encodes the tree structure. The two entries in the i-th column are the numberings
%       of the two children of the node with numbering i.
%       The root node has numbering 0, with its two children always having numbering 1 and numbering 2.
%       Thus the root node is NOT included in the 'tree' variable.
% splits: An array of length k-1. It keeps track of the numberings of the nodes being split
%         from the 1st split to the (k-1)-th split. (The first entry is always 0.)
% is_leaf: An array of length 2*(k-1). A "1" at index i means that the node with numbering i is a leaf node
%          in the final tree generated, and "0" indicates non-leaf nodes in the final tree.
% clusters: A cell array of length 2*(k-1). The i-th element contains the subset of items
%           at the node with numbering i.
% timings: An array of length k-1.
%          Its i-th element is the wall-time for performing i splits (with i+1 leaf nodes).
% Ws: A cell array of length 2*(k-1).
%     Its i-th element is the topic vector of the cluster at the node with numbering i.
% priorities: An array of length 2*(k-1).
%             Its i-th element is the modified NDCG scores at the node with numbering i (see the reference paper).
%
% Tips:
% If you want to have the flat partitioning induced by the leaf nodes in the final tree,
% use this script:
%
% partitioning = zeros(1, n); % n is the total number of data points
% leaf_level = clusters(is_leaf == 1);
% for i = 1 : length(leaf_level)
%     partitioning(leaf_level{i}) = i;
% end
% 
% (Now the entries in 'partitioning' having value 0 indicate outliers that do not belong to any cluster.)
%
%
% References:
% Da Kuang, Haesun Park
% Fast rank-2 nonnegative matrix factorization for hierarchical document clustering, KDD 2013
%
% Feb 2013
if ~exist('params', 'var')
	trial_allowance = 3;
	unbalanced = 0.1;
	vec_norm = 2.0;
	normW = true;
	anls_alg = @anls_entry_rank2_precompute;
	tol = 1e-4;
	maxiter = 10000;
else
	if isfield(params, 'trial_allowance')
		trial_allowance = params.trial_allowance;
	else
		trial_allowance = 3;
	end
	if isfield(params, 'unbalanced')
		unbalanced = params.unbalanced;
	else
		unbalanced = 0.1;
	end
	if isfield(params, 'vec_norm')
		vec_norm = params.vec_norm;
	else
		vec_norm = 2.0;
	end
	if isfield(params, 'normW')
		normW = params.normW;
	else
		normW = true;
	end
	if isfield(params, 'anls_alg')
		anls_alg = params.anls_alg;
	else
		anls_alg = @anls_entry_rank2_precompute;
	end
	if isfield(params, 'tol')
		tol = params.tol;
	else
		tol = 1e-4;
	end
	if isfield(params, 'maxiter')
		maxiter = params.maxiter;
	else
		maxiter = 10000;
	end
end

params = [];
params.vec_norm = vec_norm;
params.normW = normW;
params.anls_alg = anls_alg;
params.tol = tol;
params.maxiter = maxiter;

t0 = tic;

[m, n] = size(X);

cluster_idx=unique(cl_idx(idx));
sub_k=length(cluster_idx);
ctrary = zeros(size(X,1),sub_k);
timings = zeros(1, k-sub_k);
clusters = cell(1, 2*(k-sub_k));
Ws = cell(1, 2*(k-sub_k));
W_buffer = cell(1, 2*(k-sub_k));
H_buffer = cell(1, 2*(k-sub_k));
priorities = zeros(1, 2*(k-sub_k));
is_leaf = -1 * ones(1, 2*(k-sub_k));
tree = zeros(2, 2*(k-sub_k));
splits = -1 * ones(1, k-sub_k);
min_priority = 1e308;

% if sub_k == 1
%     term_subset = find(sum(X, 2) ~= 0);
%     W = rand(length(term_subset), 2);
%     H = rand(2, n);
%     if length(term_subset) == m
%         [W, H] = nmfsh_comb_rank2(X, W, H, params);
%     else
%         [W_tmp, H] = nmfsh_comb_rank2(X(term_subset, :), W, H, params);
%         W = zeros(m, 2);
%         W(term_subset, :) = W_tmp;
%         clear W_tmp;
%     end
% end

%initial clusters, node
for i=1:sub_k
    subidx = find(cl_idx(idx)==cluster_idx(i));
    ctrary(:,i) = mean(X(:, subidx),2);
    Ws{i} = ctrary(:,i);
    is_leaf(i) = 1;
    [~, W_buffer_one, H_buffer_one, priority_one]...
        = trial_split(trial_allowance, unbalanced, min_priority, X, subidx, Ws{i}, params);
    clusters{i} = subidx;
    W_buffer{i} = W_buffer_one;
    H_buffer{i} = H_buffer_one;
    priorities(i) = priority_one;
end

result_used = sub_k;
for i = 1 : k-sub_k
    timings(i) = toc(t0);
    
    leaves = find(is_leaf == 1);
    temp_priority = priorities(leaves);
    min_priority = min(temp_priority(temp_priority > 0));
    [max_priority, split_node] = max(temp_priority);
    if max_priority < 0
        fprintf('Cannot generate all %d leaf clusters\n', k);
        return;
    end
    split_node = leaves(split_node);
    is_leaf(split_node) = 0;
    W = W_buffer{split_node};
    H = H_buffer{split_node};
    split_subset = clusters{split_node};
    new_nodes = [result_used+1 result_used+2];
    tree(1, split_node) = new_nodes(1);
    tree(2, split_node) = new_nodes(2);
    
    result_used = result_used + 2;
	[max_val, cluster_subset] = max(H);
	clusters{new_nodes(1)} = split_subset(find(cluster_subset == 1));
	clusters{new_nodes(2)} = split_subset(find(cluster_subset == 2));
	Ws{new_nodes(1)} = W(:, 1);
	Ws{new_nodes(2)} = W(:, 2);
	splits(i) = split_node;
	is_leaf(new_nodes) = 1;

	subset = clusters{new_nodes(1)};
	[subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, X, subset, W(:, 1), params);
	clusters{new_nodes(1)} = subset;
	W_buffer{new_nodes(1)} = W_buffer_one;
	H_buffer{new_nodes(1)} = H_buffer_one;
	priorities(new_nodes(1)) = priority_one;

	subset = clusters{new_nodes(2)};
	[subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, X, subset, W(:, 2), params);
	clusters{new_nodes(2)} = subset;
	W_buffer{new_nodes(2)} = W_buffer_one;
	H_buffer{new_nodes(2)} = H_buffer_one;
	priorities(new_nodes(2)) = priority_one;
end

W = cell2mat(Ws(find(is_leaf)));
[H,temp,suc_H,numChol_H,numEq_H] = nnlsm_blockpivot(W'*W,W'*X,1,rand(k,size(X,2)));

%--------------------------------------

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

%--------------------------------------

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

%--------------------------------------

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

%--------------------------------------

function score = NDCG_part(ground, test, weight, weight_part)

[sorted, seq_idx] = sort(ground);
weight_part = weight_part(seq_idx);

n = length(test);
uncum_score = weight_part(test);
uncum_score(2:n) = uncum_score(2:n) ./ log2(2:n)';
cum_score = cumsum(uncum_score);

ideal_score = sort(weight, 'descend');
ideal_score(2:n) = ideal_score(2:n) ./ log2(2:n)';
cum_ideal_score = cumsum(ideal_score);

score = cum_score ./ cum_ideal_score;
score = score(end);
