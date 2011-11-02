function stringfile=read_stringfile(filename,whitespace,separators);
% Reads a text file into a (nlines * 1) cell array of (k * 1) cell
% arrays where k can be different for each line:
% whitespace separated substrings in the line are
% converted into a cell array of the substrings.
% Inputs: 
%   filename: name of the file to be read
%   whitespace: characters that lie between substrings. Several of these can be between two substrings.
%   separators: characters that lie between substrings. Only one such character can be between two substrings.
% The difference between 'whitespace' and 'separator' is in strings like
% "A;;B". If ";" is a whitespace character, then the example line has two
% substrings "A" and "B". If ";" is a separator character, then the 
% example line has three substrings "A", "", and "B".

% Author: Jaakko Peltonen, Jan 17, 2011 (based on old code from
% March 29, 2006) 

if nargin<3,
  separators=[];
end;

stringfile={};

f=fopen(filename,'r');
tline=fgets(f);

lines_read=0;
max_lines=0;

while tline~=-1,
  lines_read=lines_read+1;
  if mod(lines_read,1000)==0,
    lines_read
  end;
  
  a={};
  i=1; j=1;
  while i <= length(tline),
    nosubstring=1;

    % Find the beginning of a content substring (content, not whitespace)
    while nosubstring==1,
      if sum(double(tline(i))==double(whitespace))==0, 
	nosubstring=0; 
      else 
	i=i+1; 
	if i > length(tline), 
	  nosubstring=0; 
	end; 
      end;
    end;

    % If a substring was found, extract it from the line
    if i <= length(tline),
      % Find the end of the content string
      j=i;      
      while (j <= length(tline)) && (sum(double(tline(j))==double([whitespace separators]))==0), 
	j=j+1; 
      end;
      substring=tline(i:j-1);
      a={a{:},substring};
      i=j+1;
    end;
  end;

  if lines_read>max_lines,
    max_lines=lines_read*2;
    stringfile2=cell(max_lines,1);
    for k=1:(lines_read-1), stringfile2{k}=stringfile{k}; end;
    stringfile=stringfile2; stringfile2=[];
  end;

  stringfile{lines_read}=a;
  tline=fgets(f);
end;

stringfile={stringfile{1:lines_read}};

fclose(f);
