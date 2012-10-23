function f = fermidirac(ev,efermi,Tbeta);
%
% usage: f = fermidirac(ev, efermi, Tbeta); 
%
f = 1./(1+exp(Tbeta*(ev-efermi)));
