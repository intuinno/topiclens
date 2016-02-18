function [str_arr,fnames] = merge_to_single_text(dir_name,out_fname)

res = dir(dir_name);

str_arr = cell(length(res)-2,1);
fnames = cell(length(res)-2,1);

% remove the two entries, . and ..
res = res(3:end);

for i=1:(length(res))
    fid=fopen([dir_name '\\' res(i).name]);
    tmp_str= [];
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        tline(tline==13 | tline==10 | tline==26) = ' ';
        tmp_str = [tmp_str ' ' tline];
    end
    tmp_str = tmp_str(2:end);
    fclose(fid);
    str_arr{i} = tmp_str;
    fnames{i} = res(i).name;
end

if exist('out_fname','var')
    fid=fopen(out_fname,'w');
    for i=1:length(str_arr)
        fprintf(fid,'%s', str_arr{i});
        if i~=length(str_arr)
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
end