function [X1,K]=EKSM_Riccati_control(A,B,C,invR,R,m,tol,period,symm)
%
%     Extended Krylov for Riccati on   x' = A x + B u
%                                      y  = C x
%
%     A' X + X A - X B invR B' X + C' C = 0
% 
% Input
%     A           stable coeff matrix
%     B           control matrix  (tall)
%     C           observation matrix (fat)
%     invR        coef matrix of second order term
%     m           max number of iterations
%     tol         final residual tolerance (abs residual norm)
%     period      how often convergence is checked (period=1 at each iteration)
%     symm        (symm=1 A symmetric)
%
% Output
% X1 such that  X1*X1' \approx  X
% K  feedback matrix.
% Note: only feedback matrix is returned if symm=1, and the approx basis is not stored

t = cputime;

nrmb=norm(C','fro')^2;

[n,n2]=size(C');
[nb1,nb2]=size(B);

beta = norm((C)');beta2=beta^2;

pp=symamd(A);
[LA,UA,pA]=lu(A(pp,pp)','vector'); pA=pA(:);

sh=size(C',2);
s=2*sh;
cc=C(:,pp)';
rhs1=UA\(LA\cc(pA,:));  %rhs1=A'\C';
%rhs1(pA,:)=rhs1;
rhs1(pp,:)=rhs1;

[U(:,1:s),beta]=qr([C', rhs1],0); ibeta=inv(beta);
cc=U(:,1:s)'*C';
projrhs=sparse(m*s,m*s); projrhs(1:s,1:s)=cc*cc'; 
H=sparse((m+1)*s,m*s);
L=sparse((m+1)*s,m);
errore=norm(full(ibeta));
U0=0*U;
Y=[];
odds=[];
X=0;Y=0;
Bm(1:s,:)=U(:,1:s)'*B;
%xm_state(1:s,:)=U(:,1:s)'*x_state;
cc=zeros(2,2);

for j=1:m,

    jms=(j-1)*s+1;j1s=(j+1)*s;js=j*s;js1=js+1; jsh=(j-1)*s+sh;
       Up(1:n,1:sh) = A'*U(:,jms:jsh); %Up(:,1) = A1*U(:,jms); %(same space as with A)
       cc=U(pp,jsh+1:js);
       Up(:,sh+1:s) = UA\(LA\cc(pA,:));% Up(:,2) = UA\(LA\U(:,js));
       Up(pA,sh+1:s) =Up(:,sh+1:s);
       Up(pp,sh+1:s) =Up(:,sh+1:s);

%new bases block (modified gram)
    if symm
         k1=(j-1)*s+1; k2=j*s;
         coef= U'*Up;
         H(k1:k2,jms:js) = coef;
        %H(k1:k2,jms:js) = H(k1:k2,jms:js)+ coef;
         Up = Up - U*coef;
         if (j>1),
             k1=(j-2)*s+1; k2=(j-1)*s;
             coef= U0'*Up;
             H(k1:k2,jms:js) = coef;   %H(k1:k2,jms:js)+ coef;
             Up = Up - U0*coef;
         end
    else
        for il=1:2
            for kk=1:j
                k1=(kk-1)*s+1; k2=kk*s;
                coef= U(:,k1:k2)'*Up;
                H(k1:k2,jms:js) = H(k1:k2,jms:js)+ coef; % U(1:n,k1:k2)'*Up;
                Up = Up - U(:,k1:k2)*coef; %H(k1:k2,jms:js);
            end
        end
    end

% if (j<=m)
       [ww,hh]=qr(Up,0); hinv=inv(hh);
       if symm
            U0 = U;
            U  = ww;
            Bm(js1:j1s,:)=U'*B;
       %    xm_state(js1:j1s,:)=U'*x_state;
       else
            U(:,js1:j1s)=ww;
            Bm(js1:j1s,:)=U(:,js1:j1s)'*B;
       %    xm_state(js1:j1s,:)=U(:,js1:j1s)'*x_state;
       end
       H(js1:j1s,jms:js)=hh;

       I=speye(js+s);
       if (j==1),
           L(1:j*s+sh,(j-1)*sh+1:j*sh) =...
           [ H(1:s+sh,1:sh)*beta(1:sh,1:sh), speye(s+sh,sh)*beta(1:sh,1:sh)]*ibeta(1:s,sh+1:s);
       else
           L(1:j*s+s,(j-1)*sh+1:j*sh) = L(1:j*s+s,(j-1)*sh+1:j*sh) + H(1:j*s+s,jms:jms-1+sh)*rho;
       end

       odds = [odds, jms:(jms-1+sh)];   % store the odd block columns
       evens = 1:js; evens(odds)=[];
       T(1:js+s,odds)=H(1:js+s,odds);   %odd columns

       T(1:js+sh,evens)=L(1:js+sh,1:j*sh);   %even columns
       L(max(1,(j-2)*s+1):j*s+s,j*sh+1:(j+1)*sh) = ...
            ( I(max(1,(j-2)*s+1):j*s+s,(js-sh+1):js)- T(max(1,(j-2)*s+1):js+s,1:js)*H(1:js,js-sh+1:js))*hinv(sh+1:s,sh+1:s);
       rho = hinv(1:sh,1:sh)\hinv(1:sh,sh+1:s);

       I = speye(j*s);
       k=j;

if symm
   if j==1
    T(js+1:j1s,1:j1s)=Up'*A*[U,Up];
    T(1:j1s,js+1:j1s)= T(js+1:j1s,1:j1s)';
   else
    T(js+1:j1s,(j-2)*s+1:j1s)=Up'*A*[U0,U,Up];
    T((j-2)*s+1:j1s,js+1:j1s)=T(js+1:j1s,(j-2)*s+1:j1s)';
   end
end

     if (rem(j,period)==0)
       [Y,kk,ll,info] = icare(full(T(1:js,1:js))',full(Bm(1:js,:)),full(projrhs(1:js,1:js)),R); 
%keyboard
% True residual
% X = U(:,1:js)*Y*U(:,1:js)';
% norm(A'*X+X*A-X*B/R*B'*X+C'*C,'fro')/nrmb

  %er2(k)=norm(R,'fro');  %/nrmb;

if info.Report==0,
  cc = [H(js1:j1s,js-s+1:js-sh), L(js1:j1s,(j-1)*sh+1:j*sh)];
  er2(j)=sqrt(2)*norm(cc*Y(js-s+1:js,:),'fro')/nrmb;
else
  er2(j)=er2(j-1); fprintf('Warning: Reduced Riccati solution not found\n');
end

  tot_time(j)=  cputime-t;

  if (er2(j)<tol), 
     break
  end
end
end

  [uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
  sY=flipud(sY); uY=uY(:,id(end:-1:1));
  is=sum(abs(sY)>1e-12);
  Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
  Kred = + R\(Bm(1:js,:)'*Y);
% Second EKSM pass
  if symm
     [U(:,1:s),~]=qr([C', rhs1],0);
     K = Kred(:,1:s)*U(:,1:s)';
     for j1=1:j-1
        jms=(j1-1)*s+1;j1s=(j1+1)*s;js=j1*s;js1=js+1; jsh=(j1-1)*s+sh;
        Up(1:n,1:sh) = A'*U(:,1:sh); %Up(:,1) = A1*U(:,jms); %(same space as with A)
        Up(:,sh+1:s) = UA\(LA\U(pA,sh+1:s));% Up(:,2) = UA\(LA\U(:,js));
        Up(pA,sh+1:s) = Up(:,sh+1:s);
        k1=(j1-1)*s+1; k2=j1*s;
        Up = Up - U*H(k1:k2,jms:js);
        if (j1>1),
           k1=(j1-2)*s+1; k2=(j1-1)*s;
           Up = Up - U0*H(k1:k2,jms:js);
        end
        [Up,~]=qr(Up,0);
       %Up=Up/H(js1:j1s,jms:js);
       %svd(full(H(js1:j1s,jms:js)))
        K = K + Kred(:,js1:j1s)*Up';
        U0=U; U=Up;
      end
     
   else
     K = Kred*U(:,1:js)';
   end


  if symm
     %fprintf('no solution returned, only feedback matrix\n');
     X1=0;
  else
     X1 = U(:,1:js)*Y0; 
  end
space_dim = js;
rank_sol  = is;
 %er2(j+1) = er2(j); %norm(cc*Y(js-1:js,:));
  tot_time(j+1)=  cputime-t;
  total_time=cputime-t;
 nrmR = er2(end);









