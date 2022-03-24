% -------------------------------------------------------------------------
function [c,f,s] = prespde(x,t,u,DuDx,k,Pi,Qw,n,ni,Beta)
c = (ni*Beta*n)/k;
f = DuDx;
s = 0;
% -------------------------------------------------------------------------