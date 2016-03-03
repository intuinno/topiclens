load 'mnist_train.mat'
ind = randperm(size(train_X, 1));
nsample = 500;
train_X = train_X(ind(1:nsample),:);
train_labels = train_labels(ind(1:nsample));

% Set parameters
no_dims = 2;
initial_dims = 50;
perplexity = 30;
% Run t?SNE
mappedX = tsne(train_X, train_labels, no_dims, initial_dims, perplexity);

%%
% Plot results
gscatter(mappedX(:,1), mappedX(:,2), train_labels);