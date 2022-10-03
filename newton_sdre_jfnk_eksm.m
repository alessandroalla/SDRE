function [y_new,u_new] = newton_sdre_jfnk_eksm(t,ic,full_sol,PDE,g,test)


[LA,UA,pA]=lu(PDE.A','vector');pA=pA(:);
symm = 0;
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
    
      AA=PDE.A+PDE.mu2*spdiags(xold,0:0,nA,nA)-PDE.mu3*spdiags(xold.^2,0:0,nA,nA);
    
      [~,K]=EKSM_Riccati_control(AA,PDE.B,full(PDE.C_tmp),R1,full(inv(R1)),params.m,params.tol,1,symm);
    
      rhs_full = full_sol(x_sol,xnew,K,PDE,t(i),g);
      tmp_jac = @(v) jacobiano(full_sol,x_sol,xnew,K,PDE,t(i),g,v,rhs_full);
      [tmp,flag,res_gmres,iter_gmres] = gmres(tmp_jac,rhs_full,100,1e-6,1,LA,UA); 
   
      xnew = xold -tmp;
      err=norm(xnew-xold);
      xold = xnew;
      count = count+1;

    end

    fprintf('%d ...', count)
    y_new(:,i+1) = xnew;
    u_new(:,i) = -K*xnew;
end





end

