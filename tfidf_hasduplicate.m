function [X, term_subset, doc_subset, voc] = tfidf_hasduplicate(X, min_term_occurrence, min_doc_length)
% min_term_occurrence = 5
% min_doc_length = 2

stopwords = text2cell('spanish_roman.stop', '\t');
stopword_set = java.util.HashSet;
for i = 1 : length(stopwords)
	stopword_set.add(stopwords{i});
end

voc = text2cell('vocabulary.txt', '\t');

[m, n] = size(X);
term_subset = find(sum(X~=0, 2) >= min_term_occurrence & sum(X~=0, 2) < n);
doc_subset = find(sum(X) >= min_doc_length);

X = X(:, doc_subset);
%[sorted, I, J] = unique(X', 'rows');
%if length(I) ~= size(X, 2)
%	subset = sort(I);
%	X = X(:, subset);
%	doc_subset = doc_subset(subset);
%end

X = X(term_subset, :);
subset = true(1, length(term_subset));
for i = 1 : length(term_subset)
	if stopword_set.contains(voc{term_subset(i)})
		subset(i) = false;
	end
end
X = X(subset, :);
term_subset = term_subset(subset);
voc = voc(term_subset);
cell2text(voc, 'vocabulary_subset.txt');

subset1 = find(sum(X, 2) ~= 0);
subset2 = find(sum(X) ~= 0);
X = X(subset1, subset2);
term_subset = term_subset(subset1);
doc_subset = doc_subset(subset2);
voc = voc(subset1);

%pattern = (X ~= 0);
%df = full(sum(pattern, 2));
%idf = log(n./df);

%[idx, jdx, vals] = find(X);
%X = sparse(idx, jdx, log(vals) + 1);

%X = bsxfun(@times, X, idf);

D = full(1./sqrt(sum(X.^2)));
X = bsxfun(@times, X, D);
