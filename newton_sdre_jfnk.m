function [y_new,u_new,K_store] = newton_sdre_jfnk(t,ic,full_sol,PDE,g,test)

[LA,UA,pA]=lu(PDE.A','vector');pA=pA(:);
y_new(:,1) = ic(:);
params = param_Ricc;
R1=inv(PDE.R) - PDE.inf*inv(PDE.P)/(2*PDE.gamma^2);
R_gamma = inv(PDE.R) + PDE.inf*inv(PDE.P)/(PDE.gamma^2);
RR = PDE.R - PDE.inf*PDE.P/(2*PDE.gamma^2);

nA = size(PDE.A,1);
for i = 1:length(t)-1
     err = 10;
    x_sol = y_new(:,i);
    xnew = y_new(:,i);
    xold = x_sol;
    count = 0;
    while(err>1e-6)  
        switch test
            case {1,2,3}
       AA=PDE.A+PDE.mu2*spdiags(xold,0:0,nA,nA)-PDE.mu3*spdiags(xold.^2,0:0,nA,nA);                   
            case {4,5}
              [~, TT] =f_adv_new(PDE,xold);     
              AA = PDE.A - TT+1.5*spdiags(exp(-.1*xold),0:0,nA,nA);    
        end
       [Z,~,~,~,~] = ARKSM_Riccati_indef(AA, PDE.B, PDE.C_tmp,R1,params);                    
       Z=real(Z);
       K = ((((R_gamma)*PDE.B')*Z)*Z');

       rhs_full = full_sol(x_sol,xnew,K,PDE,t(i),g);
       tmp_jac = @(v) jacobiano(full_sol,x_sol,xnew,K,PDE,t(i),g,v,rhs_full);
       [tmp,flag,res_gmres,iter_gmres] = gmres(tmp_jac,rhs_full,100,1e-6,1,LA,UA); %,sol_old);
   
       xnew = xold -tmp;
       err=norm(xnew-xold);
       xold = xnew;
       count = count+1;
    end
%   K_store(:,i) = K;
    K_store=[];
    fprintf('%d ...', count)
    y_new(:,i+1) = xnew;
    u_new(:,i) = -K*xnew;
end





end

