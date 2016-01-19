function [] = writesparsetext( input,vocab,outputfilename )
%input is a sparse bag of word matrix.
%every column is a document and rows are vocab index.
[~,numDocs]=size(input);
fid=fopen(outputfilename,'wt');
%outputCell = cell(numDocs,1);
for i=1:numDocs
    idxs=find(input(:,i));    
    %for j = 1:size(idxs)
    %    oneemail{j}=vocab(idxs(j));
    %end
    for j=1:size(idxs)
        temp=vocab(idxs(j));
        fprintf(fid,'%s ',temp{1});
    end
    fprintf(fid,'\n');    
end
fclose(fid);
end