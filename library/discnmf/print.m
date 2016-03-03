function topics = print(W,U,par)
    k = par.k; dic = par.dic;
    wordmaxlen = par.wordmaxlen;
    topkwords = par.topkwords;
    [val,IX] = sort([W U],'descend'); topw = IX(1:topkwords,:);

    topics = [];
    digit = floor(log10(max([par.num1;par.num2])))+1;
    num = num2str(par.num1,['%0' num2str(digit) 'd']); num2 = num2str(par.num2,['%0' num2str(digit) 'd']); 
    len = (wordmaxlen+1+size(num,2));
    for i = 1:k
        tmp_dicA = dic(topw(:,i));
        tmp_dicB = dic(topw(:,i+k));

        for j=1:size(tmp_dicA,1)
            if length(tmp_dicA{j}) < wordmaxlen
                tmp_dicA{j} = [tmp_dicA{j} repmat(' ', 1,wordmaxlen-length(tmp_dicA{j}))];
            else
                tmp_dicA{j} = tmp_dicA{j}(1:wordmaxlen);
            end
            if length(tmp_dicB{j}) < wordmaxlen
                tmp_dicB{j} = [tmp_dicB{j} repmat(' ', 1,wordmaxlen-length(tmp_dicB{j}))];
            else
                tmp_dicB{j} = tmp_dicB{j}(1:wordmaxlen);
            end
        end
        tmp_dicA = cell2mat(tmp_dicA);tmp_dicB = cell2mat(tmp_dicB);

        topics = [topics [num(topw(:,i),:),repmat(' ',topkwords,1),tmp_dicA; repmat('-',1,len); num2(topw(:,i+k),:),repmat(' ',topkwords,1),tmp_dicB]];
    end
    
    topics
    disp(W'*U)
end