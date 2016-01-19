function [] = parsenmffile(tilesfile,nmf_dir,vocab_dir,topk)
tiles=load(tilesfile);
for i=1:numel(tiles)
    %construct file names
    hfilename=[nmf_dir '/htile' num2str(tiles(i)) '.mtx.csv'];
    wfilename=[nmf_dir '/wtile' num2str(tiles(i)) '.mtx.csv'];
    woutputfile=[nmf_dir '/wtopk' num2str(tiles(i)) '.mtx.csv'];
    houtputfile=[nmf_dir '/htopk' num2str(tiles(i)) '.mtx.csv'];
    vocabfilename=[vocab_dir '/vocabtile' num2str(tiles(i))];
    %load nmf files
    W=load(wfilename);
    H=load(hfilename);
    [m,k]=size(W);
    [~,n]=size(H);    
    %load vocabulary file
    fid = fopen(vocabfilename);
    vocab = textscan(fid,'%s','delimiter','\n');    
    vocab = vocab{1};    
    for i=1:length(vocab)        
        C = textscan(vocab{i},'%s');
        C = C{1};
        vocab{i}=C{1};
    end
    fclose(fid);
    %construct output matrices.
    Wtopk=cell(k,topk);
    Htopk=zeros(n,topk);
    %generate topk words based on W
    for i=1:k
        [~,idx]=sort(W(:,i),'descend');
        for j=1:topk
            Wtopk{i,j}=vocab{idx(j)};
        end        
    end
    %find top k topics for every document/tweet based on H
    for j=1:n
        [~,idx]=sort(H(:,j),'descend');
        Htopk(j,:)=idx(1:topk);
    end
    csvwritecellofstring(woutputfile,Wtopk);
    csvwrite(houtputfile,Htopk);
end
end