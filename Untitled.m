%%
X = rand(200,50);

[tree, splits, is_leaf, clusters, timings, Ws, priorities] = hier8_neat(X, 10);

%%

X = rand(200,50);

[tree, splits, is_leaf, clusters, timings, Ws, priorities, W, H] = hier8_neat(X, 10);
