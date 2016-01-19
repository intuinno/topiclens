function[]=csvwritecellofstring(filename,x)
fid=fopen(filename,'wt');
[rows,~]=size(x);
for i=1:rows
      fprintf(fid,'%s,',x{i,1:end-1});
      fprintf(fid,'%s\n',x{i,end});
end
fclose(fid);
end