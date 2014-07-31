function str=smallprint(mynum,maxdigits);

done=0;
if mynum==0,
  str='0';
  done=1;
end;
if (abs(mynum)>1) && (abs(mynum)<10^(maxdigits+3)) && (abs(mynum-round(mynum))<10^(-maxdigits)),
abs(mynum)
abs(mynum-round(mynum))
  str=sprintf('%d',round(mynum));
  done=1;
end;

if (done==0),
  formatstring=sprintf('%%.%de',maxdigits);
  str1=sprintf(formatstring,mynum);
  % remove useless zeroes before or after exponent
  k1=find(str1=='.');
  k2=find(str1=='e');
  str1b=str1(1:k1-1);
  uselesszeros=0;
  for l=k1+maxdigits:-1:k1+1,
    if (str1(l)=='0'),
      uselesszeros=uselesszeros+1;
    end;
  end;
  uselesszeros;
  if uselesszeros<maxdigits,
    str1b=[str1b str1(k1:k1+maxdigits-uselesszeros)];
  end;
  exponent=str2double(str1(k2+1:end));
  exponent;
  if exponent~=0,
    str1b=[str1b sprintf('e%d',exponent)];
  end;
  str=str1b;
  
  if abs(mynum)<1,
    musthavezeros=ceil(-log10(abs(mynum)))-1;
    formatstring2=sprintf('%%.%df',musthavezeros+maxdigits);
    str2=sprintf(formatstring2,mynum);
    if length(str2)<length(str),
      str=str2;
    end;
  end;
  
  done=1;
end;

