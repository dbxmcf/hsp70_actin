prologues:=3;
verbatimtex
%&latex
%\documentclass{minimal}
%\begin{document}
etex
outputtemplate := "%j.ps";
beginfig (1);
  org = (5cm,5cm);
  a=5cm;b=20cm;
  c=4;
  z0 = (a,c*a)+(5cm,5cm); z1 = (c*a,c*a)+(5cm,5cm);
  z2 = (c*a,a)+(5cm,5cm);
  draw z0--z1--z2--cycle;
  %dotlabel.top("(0,0)", z0);
  %label.bot(btex $\lbrace$ etex rotated 90
  %xscaled 5 yscaled 1.4, 1/2[z0,z1]);
  %label.rt((btex $\rbrace$ etex) xscaled 1.3
  %yscaled 3, 1/2[z1,z2]);
  %label(btex $\lbrace$ etex xscaled 1.5 yscaled 5.7
  %rotated -60, 1/2[z0,z2] + dir(120)*2mm);
  %labeloffset:=3.5mm;
  %label.bot(btex $\sqrt{3}$ etex, 1/2[z0,z1]);
  %label.rt(btex 1 etex, 1/2[z1,z2]);
  %label(btex 2 etex, 1/2[z0,z2]+dir(120)*5mm);
endfig;
end.
