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
end