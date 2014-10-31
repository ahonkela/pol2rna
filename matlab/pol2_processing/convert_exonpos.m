addpath /share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/pol2rnaseq/matlab

C=read_stringfile('transcript_exonpos.txt',[',' 10 13]);
A=zeros(length(C),length(C{2}));
for i=2:length(C),
  if mod(i,1000)==0,
    i
  end;
  A(i-1,1)=str2double(C{i}{1}(5:end));
  A(i-1,2)=str2double(C{i}{2}(5:end));
  A(i-1,3)=str2double(C{i}{3});
  A(i-1,4)=str2double(C{i}{4});
  A(i-1,5)=str2double(C{i}{5});
end;

save transcript_exonpos.mat A -mat
