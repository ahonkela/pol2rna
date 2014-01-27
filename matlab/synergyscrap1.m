;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

f=fopen('mismatch_oldnewgenes.txt','w');
for k=1:size(bininfo_new,1),
  Itemp=find(bininfo(:,5)==bininfo_new(k,5));
  if length(Itemp)==1,
    allok=1;
    if bininfo(Itemp,1)~=bininfo_new(k,1),
      allok=0;
    end;
    if bininfo(Itemp,2)~=bininfo_new(k,2),
      allok=0;
    end;
    if bininfo(Itemp,3)~=bininfo_new(k,3),
      allok=0;
    end;
    if bininfo(Itemp,6)~=bininfo_new(k,6),
      allok=0;
    end;
    if allok==0,
      fprintf(f,'ENSG%011d (%d %d) (%d %d) (%d %d) (%d %d)\n', ...
       bininfo_new(k,5), ...
       bininfo(Itemp,1), bininfo_new(k,1), ...
       bininfo(Itemp,2), bininfo_new(k,2), ...
       bininfo(Itemp,3), bininfo_new(k,3), ...
       bininfo(Itemp,6), bininfo_new(k,6) ...
      );
    end;
  end;
end;
fclose(f);

