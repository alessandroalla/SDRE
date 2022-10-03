function [X,XD,js,total_time]=accel_lyap_mrhs1_indefrhs(LA,UA,pA,PDE,LE,Zmat,rhs,Drhs,m,tol)
%function [X,XD,js]=accel_lyap_mrhs1_indefrhs(A,LE,rhs,Drhs,m,tol)
%
%
%
%%     Lyap for  E x' = A x + b   (or  M x' = N x + b)
%
%      E = LE LE'
%
%  rhs = b 
%  m = max space dimension, say sqrt(size(A))
%  tol = max final accuracy (in terms of relative residual)
%
%  X = solution factor   Z = X * XD* X'
%
% NOT FOR DISTRIBUTION. TEST VERSION
%

t=cputime;
A=PDE.A;
%spy(abs(A)>1e-10)
rhs=LE\rhs;
nrmb=norm(rhs,'fro')^2;  
nrmA=norm(A,'fro');
%singE=condest(LE*LE')/norm(LE*LE','fro');

% WORKING FOR PDE.B=PDE.H;
%Lft = (PDE.A -   ((PDE.B*inv(PDE.R))*(PDE.B'*Z))*Z' + ...
%                (PDE.gamma_eq*PDE.inf)*((PDE.H*inv(PDE.P))*(PDE.H'*Z))*Z')';

coef=( (inv(PDE.R)- PDE.gamma_eq*(PDE.inf)*inv(PDE.P))*(PDE.B'*Zmat) )';
extra=@(X)( -Zmat*(coef*(PDE.B'*X))  );
%extra=@(X)( -PDE.B*(coef*(Zmat'*X))  );

[n,sh]=size(rhs);
%rhs2=rhs*rhs';

Y=[];
odds=[];
%keyboard
% if norm(A-A',1)<1e-14, %%AA: costoso
%     disp('qui')
% %    UA = chol(-A); LA = -UA';
%      k_max =2 ;
% else
% %    [LA,UA]=lu(A);
%      k_max = m;
% end

% Always nonsym operator
      k_max = m;

Ism=speye(size(coef,1));
s=2*sh;
wrk=LE*rhs;
% Solve system (A + (-Z)*coef*B') rhs1=wrk using Sherman-Morrison formula
%  rhs1=  A\wrk - A\(-Z)*inv(I+coef*B'*(A\(-Z)) ) *coef*B'*(A\wrk)
%Lft = (PDE.A -   ((PDE.B*inv(PDE.R))*(PDE.B'*Z))*Z' + ...
%                (PDE.gamma_eq*PDE.inf)*((PDE.H*inv(PDE.P))*(PDE.H'*Z))*Z')';
rhs1(pA,:)=UA\(LA\wrk(pA,:));  % rhs1=UA\(LA\rhs);
wrkSM(pA,:)=-UA\(LA\Zmat(pA,:));  % rhs1=UA\(LA\rhs);
wrkSM1=(Ism+coef*PDE.B'*wrkSM)\(coef*PDE.B');
%rhs1=LE'*(rhs1-wrkSM*(  (Ism+coef*PDE.B'*wrkSM)\(coef*PDE.B'*rhs1)) );
rhs1=A\rhs-(-A\Zmat)*(  wrkSM1*(A\rhs)); %rhs1=LE'*(UA\(LA\(LE*rhs)));  % rhs1=UA\(LA\rhs);
UU=[rhs, rhs1];

%[ibeta,U(:,1:s)]=gram_sh(UU);
[U(:,1:s),beta]=qr(UU,0);ibeta=inv(beta);
%VV = U;
H=sparse((m+1)*s,m*s);
L=sparse((m+1)*s,m*s);
beta = beta(1:sh,1:sh); 
beta2=beta*Drhs*beta';    %beta = VV'*rhs;
Kold=0*PDE.B';

for j=1:m,

    jms=(j-1)*s+1;j1s=(j+1)*s;js=j*s;js1=js+1; jsh=(j-1)*s+sh;
% Product  Up=(A+B*coef*Z')*U
    Up(:,1:sh) = LE\(A*(LE'\U(:,jms:jsh))); 
    Up(:,1:sh) = Up(:,1:sh)+LE\extra(LE'\U(:,jms:jsh));
    wrk=LE*U(:,jsh+1:js);
%   Up(pA,sh+1:s) = LE'*(UA\(LA\wrk(pA,:)));
% Solve system (A+ B*coef*Z') rhs1=wrk using Sherman-Morrison formula
    wrk1(pA,:)=UA\(LA\wrk(pA,:));  % rhs1=UA\(LA\rhs);
    Up(:,sh+1:s)=LE'*(wrk1-wrkSM*(  wrkSM1*wrk1)) ;
   %Up(:,sh+1:s) = LE'*(UA\(LA\(LE*U(:,jsh+1:js))));
  % UU=[UU,Up];
    

%new bases block (modified gram)
    for l=1:2
        k_min=max(1,j-k_max);
        for kk=k_min:j
            k1=(kk-1)*s+1; k2=kk*s;
            coefU= U(1:n,k1:k2)'*Up;
            H(k1:k2,jms:js) = H(k1:k2,jms:js)+ coefU; % U(1:n,k1:k2)'*Up;
            Up = Up - U(:,k1:k2)*coefU; %H(k1:k2,jms:js);
        end
    end

  if (j<=m)
%      [hinv,U(1:n,js1:j1s)]=gram_sh(Up);
       [U(1:n,js1:j1s),hh]=qr(Up,0);
       H(js1:j1s,jms:js)=hh; %inv(hinv); 
       hinv=inv(hh);
  end

%{
  I=speye(js+s);
  if (j==1),
      L(1:j*s+sh,(j-1)*sh+1:j*sh) =...
      [ H(1:s+sh,1:sh)*beta(1:sh,1:sh), speye(s+sh,sh)*beta(1:sh,1:sh)]*ibeta(1:s,sh+1:s);
     %[ H(1:s+sh,1:sh)/ibeta(1:sh,1:sh), speye(s+sh,sh)/ibeta(1:sh,1:sh)]*ibeta(1:s,sh+1:s);
  else
      L(1:j*s+s,(j-1)*sh+1:j*sh) = L(1:j*s+s,(j-1)*sh+1:j*sh) + H(1:j*s+s,jms:jms-1+sh)*rho;
  end

  odds = [odds, jms:(jms-1+sh)];   % store the odd block columns
  evens = 1:js; evens(odds)=[];
  T(1:js+s,odds)=H(1:js+s,odds);   %odd columns

  T(1:js+sh,evens)=L(1:js+sh,1:j*sh);   %even columns
  L(1:j*s+s,j*sh+1:(j+1)*sh) = ...
       ( I(1:j*s+s,(js-sh+1):js)- T(1:js+s,1:js)*H(1:js,js-sh+1:js))*hinv(sh+1:s,sh+1:s);
  rho = hinv(1:sh,1:sh)\hinv(1:sh,sh+1:s);

  I = speye(j*s); 
 %VV=U(:,1:j*s);
  k=j;
%}

  k=j;
% This could be made cheaper...
 T(1:js,1:js)=U(:,1:js)'*(LE\(A*(LE'\U(:,1:js))));
    T(1:js,1:js) = T(1:js,1:js)+U(:,1:js)'*(LE\extra(LE'\U(:,1:js)));
%T(1:js,1:js)=U(:,1:js)'*A*U(:,1:js);

  Y = lyap(full(T(1:js,1:js)),speye(js,sh)*beta2*speye(js,sh)');
 %Y = lyap(full(T(1:js,1:js)),speye(k*s,sh)*beta2*speye(k*s,sh)');
% nrmx=norm(Y,'fro');

  cc=full(U(:,js1:j1s)'*( LE\(A*(LE'\U(:,js-s+1:js))) + LE\extra(LE'\U(:,js-s+1:js))) );
% cc = [H(js1:j1s,js-s+1:js-sh), L(js1:j1s,(j-1)*sh+1:j*sh)];
% er2(j)=norm(cc*Y(js-s+1:js,:),'fro');  %/(nrmb+nrmA*singE*nrmx);

   er2(j) = norm(cc*Y(js-s+1:js,:),'fro'); %/(nrmb+nrma*nrmx);
 %K=(PDE.B'*U(:,1:js)*Y)*U(:,1:js)';
 % disp([k,er2(j),js,norm(K-Kold,'fro')/norm(PDE.B),norm(K-Kold,'fro')/norm(K,'fro'),norm(K,'fro')])
 %disp([j,er2(j),js])
  %Kold=K;

 % X=U(:,1:js)*Y*U(:,1:js)';
 %op1=extra(X);
 %op2=extra(X');
 % normlyap=norm( A'*X+op1 +X*A+op2' + rhs*Drhs*rhs','fro')
%er2(k)=normlyap;
  tot_time(j)=  cputime-t;
% disp([k,er2(k),js])

  if (er2(j)<tol), 
  break
  end
end
 % disp([' '])
 %  disp([j,er2(j),js]);

  [uY,sY]=eig(Y); 
%[ssY,id]=sort(diag(sY));
%  ssY=flipud(ssY); uY=uY(:,id(end:-1:1));
  ssY=diag(sY);
  trunc_tol=1e-10;
  is=find(abs(ssY)>trunc_tol);
  Y0 = uY(:,is); %*diag(sqrt(sY(1:is)));
 %Y0 = uY(:,1:is)*diag(sqrt(sY(1:is)));
  X = U(:,1:js)*Y0; XD=diag(ssY(is));
  er2(j+1) = er2(j);
  tot_time(j+1)=  cputime-t;
  total_time=cputime-t;

return
