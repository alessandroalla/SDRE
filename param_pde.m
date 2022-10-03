function [PDE,g,t_grid,dt,ic,full_sol,x] = param_pde(n,test)

% provide the parameters for the control problem 
% and the PDE 

%       State equation (SE):
%              y_t = a Δy + mu_1 y + mu_2 yˆ2- mu_3 yˆ3 + Bu + Hg
%             plus Boundary conditions and Initial condition
%
%      Running cost: ||Cy||ˆ2 + R ||u||ˆ2
%               
% Input:
%    n   # space FD discretization nodes 
% 
% Output:
%     PDE.n        # space FD discretization nodes
%     PDE.mu1      parameter mu_1 in (SE)
%     PDE.mu2      parameter mu_2 in (SE)
%     PDE.mu3      parameter mu_3 in (SE)
%     PDE.A        discretization of the operator Δ in (SE)
%     PDE.B        vector B in (SE)
%     PDE.R        paramter R in the running cost
%     PDE.H        vector H in (SE)
%     PDE.gamma    parameter gamma such that 
%                     ||Cy||ˆ2 + R ||u||ˆ2 <= gammaˆ2 ||g||ˆ2
%     PDE.inf      boolean to activate the Hinf problem
%     PDE.cal      boolean to activate H2 control within the Hinf problem
%     g            disturbance of the controlled problem
%     t_grid       discrete temporal grid
%     dt           temporal step size
%     full_sol     implicit Euler scheme for (SE)
%     x            discrete spatial domain

PDE.n = n;
dt = 0.05;
tend = 3;
t_grid = 0:dt:tend;
x = linspace(0,1,n);
y = x;
ic = sin(pi*x); % initial condition
ic = ic'*ic;


switch test
    case {1,2}
        PDE.mu1 = -0.1;
        PDE.mu2 = 10; 
        PDE.mu3 = 10; 
        a =0.2; % diffusion coeff
        [~,~,PDE.A] = laplacian([n,n],{'NN','NN'});
        PDE.A=-(a*(n-1)^2)*PDE.A+PDE.mu1*speye(size(PDE.A));
    case 3
PDE.mu1 = 0;
PDE.mu2 = 8; 
PDE.mu3 = 8; 
a =0.1; % diffusion coeff
A0=spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);
PDE.A=kron(A0,speye(n))+kron(speye(n),A0);
PDE.A=-(a*(n-1)^2)*PDE.A+PDE.mu1*speye(size(PDE.A));
    case {4,5}
%PDE.mu1 = 0;
%PDE.mu2 = 0; 
%PDE.mu3 = 0; 
a =0.1; % diffusion coeff
A0=spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);
PDE.A=kron(A0,speye(n))+kron(speye(n),A0);
PDE.A=-(a*(n-1)^2)*PDE.A;
[PDE.Tpos,PDE.Tneg] = get_adv(n); 
        
end
 
 switch test
     case {1,2,3}
  full_sol =@(y,tmp,K,PDE,t,g)(tmp-y-dt*(PDE.A*tmp+PDE.mu2*tmp.^2-PDE.mu3*tmp.^3-PDE.B*(K*tmp)-PDE.cal*PDE.H*g(t)));
     case {4,5}
  full_sol =@(y,tmp,K,PDE,t,g)(tmp-y-dt*(PDE.A*tmp-f_adv_new(PDE,tmp)-PDE.B*(K*tmp)-PDE.cal*PDE.H*g(t)+1.5*tmp.*exp(-.1*tmp)));
 end

 
% Settings for the control B, C, P, H, gamma
[X, Y] = meshgrid(x,y);
[PDE.B,PDE.C_tmp] = get_B_C(X,Y,n);
PDE.B = sum(PDE.B,2);
dim_c = size(PDE.B,2);
PDE.R = .05*speye(dim_c);
PDE.P = 1*speye(dim_c);
PDE.H = PDE.B;
PDE.gamma = 0.5;
PDE.gamma_eq = 1/PDE.gamma^2;

% H-inf control on
switch test
    case {1,3,4}
        PDE.inf = 0;
    case {2,5}
        PDE.inf = 1;
end
PDE.cal = 1;

% Define Disturbances
kk = 0.1;
g = @(t) kk*(sin(2*t));


end

