prologues:=3;
verbatimtex
%&latex
%\documentclass{minimal}
%\begin{document}
etex
outputtemplate := "%j.ps";
beginfig (1);
  pair org;
  org = (0cm,10cm);
  a=6cm;b=20cm;
  c=3;
  z0 = (a,c*a) + org;
  z1 = (c*a,c*a) + org;
  z2 = (c*a,a) + org;
  draw z0--z1--z2--cycle;
  chunks = 4;
  for i=1 upto chunks:
    %drawarrow p shifted -(i*A/10) withcolor red;
    %show i/3;
    numeric fc;
    fc=i/chunks;
    %draw 1./c*i[z0,z2] -- 1./c*i[z1,z2];
    draw fc[z0,z2] -- fc[z1,z2];
    draw fc[z0,z1] -- fc[z0,z2];
  endfor;
  %draw 1/3[z0,z2] -- 1/3[z1,z2];
  %dotlabel.top("(0,0)", z0);
  %label.bot(btex $\lbrace$ etex rotated 90
  %xscaled 5 yscaled 1.4, 1/2[z0,z1]);
  %label.rt((btex $\rbrace$ etex) xscaled 1.3
  %yscaled 3, 1/2[z1,z2]);
  %label(btex $\lbrace$ etex xscaled 1.5 yscaled 5.7
  %rotated -60, 1/2[z0,z2] + dir(120)*2mm);
  labeloffset:=3.5mm;
  label.top(btex dataN etex, 1/2[z0,z1]);
  %label.rt((btex $\rbrace$ etex) xscaled 1.3 yscaled 3, 1/2[z1,z2]);
  label.rt((btex $\rbrace$ etex) xscaled 12 yscaled 30, 1/2[z1,z2]);
  %label.rt(btex 1 etex, 1/2[z1,z2]);
  %label(btex 2 etex, 1/2[z0,z2]+dir(120)*5mm);
endfig;
end.
