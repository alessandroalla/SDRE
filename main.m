clear
clc
close all
warning('off')


% The code aims to control the following PDEs:
%    (SE-1)       y_t = a Δy + mu_1 y + mu_2 yˆ2- mu_3 yˆ3 + Bu + Hg
%                   plus Boundary conditions and Initial condition
%
%    (SE-2)       y_t = a Δy - y∂y + 1.5y .*exp(-.1y) + Bu + Hg 
%                   plus Boundary conditions and Initial condition
%
%
%
% using three approaches based on State Dependent Riccati Equation (SDRE)

% (SE) after discretization leads to the following system
%
%           y' = A(y) y + B u + Hg
%
% A(y) is taken as the discretization of 
%               aΔy + mu_1 y + mu_2 yˆ2- mu_3 yˆ3 in (SE-1)
%
%               a Δy - \nabla y + 1.5.*exp(-.1y) in (SE-2)
%
%
% The running cost is quadratic: ||Cy||ˆ2 + R ||u||ˆ2 in the H2-case
% whereas in the Hinf-control the running cost has to satisfy
%       ||Cy||ˆ2 + R ||u||ˆ2 <= gammaˆ2 P ||g||
%
% We provide three different algorithms to approximate the SDRE problem
%
% 1) At each time iteration we solve the SDRE
%               A(y)' X(y)  +  X(y) A(y) - X(y) S X(y) + C'C=0
%    where S = B inv(R) B' - inv(2 gammaˆ2) H inv(P) H
%    and set the control u = - inv(R) B'X(y) y
%
% 2) Similarly to item 1, but we do not compute explicitely X(y) at each
% iteration. Instead, we compute directly K(y) = - inv(R) B'X(y) y
% and set the control u = -K(y) y.
%
% 3) We suppose that A(y) = A_0 + Σ A_j f_j(y)
%    and X(y) \approx X_0 + Σ X_j f_j(y)
%    where X_0 solves  A_0' X_0  +  X_0 A_0 - X_0 S X_0 + C'C=0
%    and at each iteration we solve the Lyapunov equation
%      W(y)C_0 + C_0' W(y) + Σ Q_j f_j(y) = 0
%    with W(y) = Σ X_j f_j(y), C_0 = A_0 - S X_0, Q_j = X_0 A_j + A_j X_0.
% 
% We set the control u = - inv(R)B' (X_0 + W(y)) y.
%
%
%
%
%
% Description of the variables
% 
%     test        choice of the test
%     n           dimension of the proble
%     PDE         structure containg the setting of our problem 
%                  (see param_pde.m for details of the structure)
%     g           disturbance of the controlled problem
%     t_grid      discrete temporal grid
%     dt          temporal step size
%     full_sol    implicit Euler scheme for (SE)
%     x           discrete spatial domain
% 
%     y_lqr                 Controlled state using linearized lqr
%     u_lqr                 Control computed linearized lqr
%     y_sdre_mpc            Controlled state using approach 1)
%     u_sdre_mpc            Control computed using approach 1)
%     y_sdre_mpc_eksm       Controlled state using approach 2)
%     u_sdre_mpc_eksm       Control computed using approach 2)
%     y_sdre_lyap           Controlled state using approach 3)
%     u_sdre_lyap           Control computed using approach 3)


% if you use this code (corresponding to Test 3 in the ref below) please cite:
% Alessandro Alla, Dante Kalise and Valeria Simoncini 
% State-dependent Riccati equation feedback stabilization for nonlinear PDEs
% Technical Report, June 2021 https://arxiv.org/abs/2106.07163
%
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
%FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
%COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
%IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
%CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%




% PDE settings

n =101; % Dimension of the discretized problem for each dimension
test = 5; 
            %1 for H2 with (SE-1), 
            %2 for Hinf for (SE-1), 
            %3 for H2 for (SE-1) and direct computation of the feedback gain
            %4 for H2 with (SE-2), 
            %5 for Hinf for (SE-2), 
[PDE,g,t_grid,dt,ic,full_sol,x] = param_pde(n,test); % PDE settings and control problem


%% LQR
 
disp('LQR')

% computing the ARE for the linearized problem
    R1=inv(PDE.R);
    params = param_Ricc;
   [Z,~,~,~,~] = ARKSM_Riccati_indef(PDE.A, PDE.B, PDE.C_tmp,R1,params);  
    Z=real(Z);
    K_tmp = (((R1*PDE.B')*Z)*Z'); %feedback gain matrix

    [y_lqr] = newton_jfnk(t_grid,ic,full_sol,PDE,g,K_tmp); 
    fprintf('\n')
    u_lqr = - K_tmp*y_lqr; %LQR control


cost_lqr = calc_cost(PDE,y_lqr); %evaluation of the cost functional

figure
 surf(x,x,reshape(y_lqr(:,end),n,n))
  xlabel('\xi_1')
 ylabel('\xi_2')
 set(gca,'FontSize',16)



%% SDRE-MPC
disp('SDRE-MPC Alg 3.1') % Number of the Alg referes to the arxiv report

%This functions solves the SDRE using approach 1)
tic
[y_sdre_mpc, u_sdre_mpc,~] = newton_sdre_jfnk(t_grid,ic,full_sol,PDE,g,test);
fprintf('\n')
toc
fprintf('\n')
%This function solves the SDRE using approach 2)
if test ==3
    [y_sdre_mpc_eksm, u_sdre_mpc_eksm] = newton_sdre_jfnk_eksm(t_grid,ic,full_sol,PDE,g,test);
    fprintf('\n')
end
figure
 surf(x,x,reshape(y_sdre_mpc(:,end),n,n))
 xlabel('\xi_1')
 ylabel('\xi_2')
 set(gca,'FontSize',16)
 
 cost_sdre = calc_cost(PDE,y_sdre_mpc); %evaluation of the cost functional
  
%%
 disp('SDRE-LYAP Alg 3.3') % Number of the Alg referes to the arxiv report
 
 % Reordering for test 1,2,3
 switch test
     case {1,2,3}
   pperm=symamd(PDE.A);
   PDE.A=PDE.A(pperm,pperm); 
   ic = ic(:);
   ic = ic(pperm);
   PDE.B = PDE.B(pperm);
   PDE.C_tmp = PDE.C_tmp(:,pperm);
 end

 %This function solves the SDRE using approach 3)
 tic
[y_sdre_lyap, u_sdre_lyap,~] = newton_lyap_jfnk(t_grid,ic,full_sol,PDE,g,test);
fprintf('\n'); toc;
fprintf('\n')


cost_sdre_lyap = calc_cost(PDE,y_sdre_lyap);

switch test
    case {1,2,3}
         P_pperm = sparse(length(ic),length(ic));
       for i =1:size(PDE.A,1)
           P_pperm(pperm(i),i) = 1;
       end
   
end

 figure
switch test
    case {1,2,3}
 surf(x,x,reshape((P_pperm*(y_sdre_lyap(:,end))),n,n))
    case {4,5}
  surf(x,x,reshape(((y_sdre_lyap(:,end))),n,n))
end
 xlabel('\xi_1')
 ylabel('\xi_2')
 set(gca,'FontSize',16)                               
 
 fprintf('LQR cost, %f \n', dt*cost_lqr)
 fprintf('SDRE-MPC Alg 3.1 cost, %f \n', dt*cost_sdre)
 fprintf('SDRE-LYAP Alg 3.3 cost, %f \n', dt*cost_sdre_lyap)


figure
plot(t_grid(1:end-1),u_lqr(1:end-1),t_grid(1:end-1),u_sdre_mpc,t_grid(1:end-1),u_sdre_lyap,'LineWidth',2)
title('Control Input')
grid
legend('LQR','SDRE','LYAP','Location','southeast')
set(gca,'FontSize',16)
