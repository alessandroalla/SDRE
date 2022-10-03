function [y_new,u_new,K_store] = newton_lyap_jfnk(t,ic,full_sol,PDE,g,test)


[LA,UA,pA]=lu(PDE.A','vector');pA=pA(:);

y_new(:,1) = ic(:);
params = param_Ricc;
R1=inv(PDE.R) - PDE.inf*inv(PDE.P)/(2*PDE.gamma^2);
R_gamma = inv(PDE.R) + PDE.inf*inv(PDE.P)/(PDE.gamma^2);
RR = PDE.R - PDE.inf*PDE.P/(2*PDE.gamma^2);

m_lyap=100;
tol_lyap=1e-6;


[Z,~,~,~,~] = ARKSM_Riccati_indef(PDE.A, PDE.B, PDE.C_tmp,R1,params); 

Z = real(Z);
ZZ = Z;
K0=((R_gamma*PDE.B')*Z)*Z';
nz=size(ZZ,2);
II = speye(size(PDE.A));   

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
         S1 = calc_S1(xold,PDE);   % computation of Σ A_j f_j(y)        
         case {4,5}
               [~, TT] =f_adv_new(PDE,xold);
                S1 =  -TT+1.5*spdiags(exp(-.1*xold),0:0,nA,nA);  
     end
         [QU,SU,VU]=svd([ZZ, S1*ZZ],0); nz1=sum(diag(SU)>1e-7); 
         EU=SU(1:nz1,:)*(VU'*[sparse(nz,nz),speye(nz,nz);speye(nz,nz),sparse(nz,nz)]*VU)*SU(1:nz1,:)';
         [EUv,EUl]=eig(EU); [eeu,iu]=sort(abs(diag(EUl)),'descend'); EUv=EUv(:,iu); EUl=EUl(iu,iu);
         neu=sum(eeu>1e-7); QU=QU(:,1:nz1)*EUv(:,1:neu); EU=EUl(1:neu,1:neu);

         [Z_lyap_f,ZD_lyap_f,~,~]=accel_lyap_mrhs1_indefrhs(LA,UA,pA,PDE,II,Z,QU,EU,m_lyap,tol_lyap);

         K=K0 + ((R_gamma*PDE.B'*Z_lyap_f)*ZD_lyap_f)*Z_lyap_f';
       
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
    u_new(:,i) = -K*xnew;
    %K_store(:,i) = K;
    K_store = [];
end



end

