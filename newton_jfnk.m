function [y_new] = newton_jfnk(t,ic,full_sol,PDE,g,K)
y_new(:,1) = ic(:);
[LA,UA,pA]=lu(PDE.A','vector');pA=pA(:);
for i = 1:length(t)-1
    err = 10;
    x_sol = y_new(:,i);
    xnew = y_new(:,i);
    xold = x_sol;
    count = 0;
    while(err>1e-6)        
             rhs_full = full_sol(x_sol,xnew,K,PDE,t(i),g);
             tmp_jac = @(v) jacobiano(full_sol,x_sol,xnew,K,PDE,t(i),g,v,rhs_full);
             [tmp,flag,res_gmres,iter_gmres] = gmres(tmp_jac,rhs_full,100,1e-6,1,LA,UA); %,sol_old);

        
         xnew = xold -tmp;
        err=norm(xnew-xold);
        xold = xnew;
        count = count+1;
    end
    fprintf('%d ...', count)
    y_new(:,i+1) = xnew;
end


 

end

