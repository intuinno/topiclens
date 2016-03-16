leaves = find(is_leaf_sub == 1);
temp_priority = priorities_sub(leaves);
min_priority = min(temp_priority(temp_priority > 0));
[max_priority, split_node] = max(temp_priority);

if max_priority < 0
    fprintf('Cannot generate all %d leaf clusters_sub\n', k);
    return;
end
split_node = leaves(split_node);
is_leaf_sub(split_node) = 0;
W_sub = W_buffer{split_node};
H_sub = H_buffer{split_node};
split_subset = clusters_sub{split_node};
new_nodes = [result_used+1 result_used+2];
tree_sub(1, split_node) = new_nodes(1);
tree_sub(2, split_node) = new_nodes(2);

result_used = result_used + 2;
[max_val, cluster_subset] = max(H_sub);
clusters_sub{new_nodes(1)} = split_subset(cluster_subset == 1);
clusters_sub{new_nodes(2)} = split_subset(cluster_subset == 2);
Ws_sub{new_nodes(1)} = W_sub(:, 1);
Ws_sub{new_nodes(2)} = W_sub(:, 2);
splits_sub(i) = split_node;
is_leaf_sub(new_nodes) = 1;

subset = clusters_sub{new_nodes(1)};
[subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, target_A_sub, subset, W_sub(:, 1), params);
clusters_sub{new_nodes(1)} = subset;
W_buffer{new_nodes(1)} = W_buffer_one;
H_buffer{new_nodes(1)} = H_buffer_one;
priorities_sub(new_nodes(1)) = priority_one;

subset = clusters_sub{new_nodes(2)};
[subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, target_A_sub, subset, W_sub(:, 2), params);
clusters_sub{new_nodes(2)} = subset;
W_buffer{new_nodes(2)} = W_buffer_one;
H_buffer{new_nodes(2)} = H_buffer_one;
priorities_sub(new_nodes(2)) = priority_one;

W_sub = cell2mat(Ws_sub(find(is_leaf_sub)));
[H_sub,temp,suc_H,numChol_H,numEq_H] = nnlsm_blockpivot(W_sub'*W_sub,W_sub'*target_A_sub,1,rand(i+sub_k,size(target_A_sub,2)));



%%%%%%%%%% after one iteration %%%%%%%%%
[Wtopk_sub,Htopk_sub,DocTopk_sub,Wtopk_idx_sub] = parsenmf(W_sub,H_sub,dict,topk_sub);
size(Wtopk_idx_sub)

for i= 1:3
    for j=1: size(Wtopk_idx_sub,2)
        tmp = Wtopk_idx_sub(i,j);
        if (sum(sum(tmp == Wtopk_idx_sub(1:4,:))) >= 3)
            Wtopk_idx_sub(i,j) = Wtopk_idx_sub(i+1,j);
            Wtopk_idx_sub(i+1,j) = tmp;
        end
    end
end

Wtopk_idx_sub = Wtopk_idx_sub';

[~,cl_idx_sub] = max(H_sub);











