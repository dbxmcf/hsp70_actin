prologues:=3;
verbatimtex
%&latex
\documentclass{minimal}
\begin{document}
etex
outputtemplate := "%j.ps";
beginfig (1);
  % draw a line
  fill for p=(20,30),(60,30),(300,540),(40,520):
    p-- endfor cycle withcolor 0.8white;
    for p=(20,30),(60,30),(300,540),(40,520):
        dotlabel.bot("(" & decimal(xpart(p)) &
        "," & decimal(ypart(p)) & ")", p);
        endfor;
endfig;
end.
