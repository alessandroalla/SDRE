function [params] = param_Ricc

params.m=100;
params.tol=1e-6;
%
% smin and smax should be set to an  approx to the 
% extreme eigs of the Riccati matrix -A  (with A<0)
%
params.smin=5e-1; params.smax=5e1;
params.ch=0;
params.period=1;
params.Hritz=1; 


end

