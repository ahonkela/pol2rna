function add_summitcolors_to_peaks(summitfile,originalpeakfile,outputpeakfile);

f1 = fopen(summitfile,'r');
f2 = fopen(originalpeakfile,'r');
f3 = fopen(outputpeakfile,'w');

maxsummitheight=1000;
minsummitheight=1;

tline = fgets(f1);
nlines=0;
while tline ~= -1,
  nlines=nlines+1;
  if mod(nlines,1000)==0, tline, end;

  %-------------------------------
  % Find the summit height, which is the 5h column of the summit .bed file
  % Columns are separated by whitespace or tabs
  %-------------------------------
  i = 1; 
  while (double(tline(i)) ~= double(' ')) && (double(tline(i)) ~= 9), i=i+1; end; % end of column 1
  while (double(tline(i)) == double(' ')) || (double(tline(i)) == 9), i=i+1; end; % start of column 2
  while (double(tline(i)) ~= double(' ')) && (double(tline(i)) ~= 9), i=i+1; end; % end of column 2
  while (double(tline(i)) == double(' ')) || (double(tline(i)) == 9), i=i+1; end; % start of column 3
  while (double(tline(i)) ~= double(' ')) && (double(tline(i)) ~= 9), i=i+1; end; % end of column 3
  while (double(tline(i)) == double(' ')) || (double(tline(i)) == 9), i=i+1; end; % start of column 4
  while (double(tline(i)) ~= double(' ')) && (double(tline(i)) ~= 9), i=i+1; end; % end of column 4
  while (double(tline(i)) == double(' ')) || (double(tline(i)) == 9), i=i+1; end; % start of column 5
  %tline(i:end)
  summitheight = str2double(tline(i:end));

  % truncate summit height to desired minimum and maximum
  if (summitheight<minsummitheight), summitheight = minsummitheight; end;
  if (summitheight > maxsummitheight), summitheight = maxsummitheight; end;
  propheight=(summitheight-minsummitheight)/(maxsummitheight-minsummitheight);

  %[summitheight propheight]
  %pause

  % use a logarithmic scale for the coloring
  peakcolor=255*log(propheight+1)/log(2);

  %-------------------------------
  % add the summit height and other dummy fields to the peak file
  %-------------------------------
  peakline = fgets(f2);
  % first remove the original linebreak
  for i=1:length(peakline), 
    if ((peakline(i)==10) || (peakline(i)==13)), 
      peakline(i)=' '; 
    end;
  end;
  % print the additional information to the end of the file
  fprintf(f3,'%s + 0 0 %d,%d,%d\n',peakline,0,floor(peakcolor),255-floor(peakcolor));

  tline = fgets(f1);
end; 

fclose(f1);
fclose(f2);
fclose(f3);


