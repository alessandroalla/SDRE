function [S1] = calc_S1(yold,PDE)
%
%   INPUT
%       yold   current state
%       PDE    structure of the problem
%
%  OUTPUT
%       S1     nonlinear terms in the State equation to obtain
%                A(y) = A_0 + Σ A_j f_j(y) with S1 = Σ A_j f_j(y) 
%
        tmp = (PDE.mu2*yold) -(PDE.mu3*yold.^2);        
        nt=length(tmp);
        S1 = spdiags(tmp,0,nt,nt);

end

