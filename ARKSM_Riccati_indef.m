function [Z,nrmrestot,VV,K,s,RKStotal_time]=ARKSM_Riccati_indef(A,B,C,R1,params)
%function [Z,nrmrestot,VV,K,s]=ARKSM_Riccati_indef(A,B,C,R1,params);
%      
% Approximately Solve  
%                A' X  +  X A - XB R1 B'X + C'C=0
%
% by the Rational Krylov subspace method 
% (Galerkin condition onto the Rational Krylov subspace)
%
% Input:
%
% A               coeff. matrix.  n x n
% B               second order term   n x s
% C               rhs factor   p x n
% R1              inner factor. can be indefinite
% params.m        max space dimension allowed
% params.tol      stopping tolerance 
% params.s1,params.smax estimates for real spectral interval
%                 associated with field of values of A'
%                 e.g., smax=norm(A,1); s1=smax/condest(A);
% params.ch       ch=1  complex poles  ch=0 real poles
% params.period   how often check convergence (period=1 means at each iteration)
% params.Hritz    Hritz=1 uses adaptive spectral region, Hritz=0 uses spectral
%                 region only based on projected A

%
% Output:
%
% Z    factor of approximate solution  X = Z Z'
% nrmrestot  residual norm history
% VV   orthonormal basis for the approx space
% K    VV'*A'* VV
% s    sequence of generated poles
%
% Hints:
% 2) Provide "comfortable" (loose bounds) estimates s1, emax


   m=params.m;
   tol=params.tol;
   s1=params.smin;
   emax=params.smax;
   ch=params.ch;
   period=params.period;
   Hritz=params.Hritz;


ttime=cputime;
[n,n]=size(A);
C=full(C);
p=size(C',2);
I=speye(p);O=0*I;
In=speye(n);
Lres=C';
%keyboard
%[rr,V]=gram_sh(Lres); nrmb=norm(inv(rr),'fro')^2; 
[V,rr]=qr(Lres,0); nrmb=norm(rr,'fro')^2; 
beta=V'*Lres; beta2=beta*beta';
errtot=[];

VV=V;

H=sparse(p*(m+2),p*(m+1));
nrmrestot=[];
%nrma=norm(A,'fro');

if (norm(A-A',1)<1e-14), symm=1; else symm=0;end
 newAv=A'*V;
 K=full(V'*newAv);
 s=s1(1);
 eH=eig(K);
 eHpoints = sort([s1(:)',emax]);
 snew=newpolei(eHpoints,eH,s1(1)*ones(p,1));
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 s=[s,snew];

% additional steps
cmplxflag=0;
itsinner=0;

%fprintf('   space dim   residual norm \n')
i=0;

while i < m

  i=i+1;

  paired=0;
  itp=1;
  while (paired==0),

    i1=i+1; it = 0; t=0.;
    w=V;
    wrk = (A'-snew*In)\w; 
    wrk= wrk;  
% Gram-Schmidt step
    jms=(i-1)*p+1;j1s=(i+1)*p;js=i*p;js1=js+1;
    Bm(jms:js,:)= V'*B;

    for it=1:2,
      for kk=1:i
        %keyboard
        k1=(kk-1)*p+1; k2=kk*p; 
        ww=VV(1:n,k1:k2);
        gamma=ww'*wrk;
        H(k1:k2,jms:js) = H(k1:k2,jms:js)+ gamma;
        %gamma = 0; 
        wrk = wrk - ww*gamma;
      end
    end
%keyboard
    [V,hinv]=qr(wrk,0); H(js1:j1s,jms:js)=hinv; hinv = inv(hinv);
    if (cmplxflag), snew=conj(snew); s=[s,snew];cmplxflag=0;
    newAv=A'*V;
    %D = kron(spdiag(s(2:end)),I); A:modification
    D = kron(diag(s(2:end)),I);
    g = VV'*newAv;
    g1 = g; 
    g2 = V'*A'*VV;
    g3 = V'*A'*V;
    K = [K g1; g2, g3];
    VV=[VV,V];
    i=i+1; itp=itp+1;
    else, paired=1; end
  end


    ih1=i1; ih=i;
    newAv=A'*V;
    %keyboard
    %D = kron(spdiag(s(2:end)),I);
    D = kron(diag(s(2:end)),I); % A: modification
    g = VV'*newAv;

   if (symm), K=(K+K')/2;, end

if (rem(i,period)==0)

% Solve the projected problem
     rhs2=speye(ih*p,p)*beta2*speye(ih*p,p)';

%[XX,EE]=eig(full([K',  -Bm*R1*Bm'; -rhs2, -K]));  nK=size(K,1);
%[ee,ii]=sort(real(diag(EE))); xxsort=(XX(:,ii)); 
%Y=xxsort(nK+1:end,1:nK)/xxsort(1:nK,1:nK); 
%Y=real(Y);
%keyboard

  Y = icare(full(K'),full(Bm),full(rhs2),full(inv(R1)));

if isempty(Y), Y=ones(size(K));end

 %   Y = care(full(K'),full(Bm),full(rhs2),full(inv(R1)));
    %Y = care(K',Bm,rhs2);
%[eig(Y),sort(eig(K'-Y*Bm*Bm'))]
%pause
%fv(K'-Y*Bm*Bm')
%pause

% computed residual   (exact, in exact arithmetic) cheaper computation possible
     u1=newAv-VV*g;   
     d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
     U=[-V*s(end),  d u1 ];
   % keyboard
     rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:)); % A:change here size (:,2)
% abs residual
    %nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro');
% rel residual
     nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/nrmb;
%backward error
    %nrmx = norm(Y,'fro');
    %nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+nrma*nrmx);

     nrmrestot=[nrmrestot,nrmres];
     
    
     

%      disp([i,nrmres]) %AA: commenta
     if (nrmres<tol), break,end
end


% New poles and zeros
  if Hritz
    eH=sort(eig(K'-Y*Bm*R1*Bm'));
    [ih,hh]=find(real(eH)>0);
    eH(ih)=-abs(real(eH(ih)))+i*imag(eH(ih));
  else
    eH=sort(eig(K));
  end
  eHorig=eH;
  eK=sort(eig(K));

 if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH)) ~=0 & max(abs(imag(eH)))>1e-5 & length(eH)>2) % Roots lambdas come from convex hull too
     eH=[eH;-emax];
      ij=convhull(real(eH),imag(eH)); eH=eH(ij);
      ieH=length(eH); missing=ih*p-ieH;
      while missing>0,                         % include enough points from the border
        neweH=(eH(1:ieH-1)+eH(2:ieH))/2;missing=ih*p-length(eH);
        eH=[eH;neweH];
      end
    % eH=eH(1:ih);
      eHpoints=-eH;
      eH=eHorig;
    else                                  % if all real eigs, no convex hull possible
      eHpoints = sort([s1; emax.';-real(eH)]);
    end


 else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     if (any(imag(eH)) ~=0 & length(eH)>2)    % Roots lambdas come from convex hull too
       eH=[eH;-s1;-emax.'];
       ij=convhull(real(eH),imag(eH)); eH=eH(ij);
       ieH=length(eH); missing=ih*p-ieH;
       while missing>0, % include enough points from the border
         neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
         eH=[eH;neweH];
         missing=ih*p-length(eH);
       end
       eH=eH(1:ih*p);
     end
      eHpoints = sort([s1; emax.';-real(eH)]);
      eH=eHorig;
 end


 gs=kron(s(2:end),ones(1,p))';

 snew = newpolei(eHpoints,eH,gs);
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end


% If pole is complex, include its conjugate
 if (imag(snew) ~=0), cmplxflag=1;end
 s=[s,snew];

 g1 = g; 
 g2 = V'*A'*VV;
 g3 = V'*A'*V;
 K = [K g1; g2, g3];
 VV=[VV,V];

end;
      %disp([i,nrmres])

 %disp([i,nrmres]) %AA: commenta

% factored solution 
[uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
sY=flipud(sY); uY=uY(:,id(end:-1:1));
is=sum(abs(sY)>1e-15);
Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
Z = VV(:,1:size(Y0,1))*Y0; 

final_rank=is;
RKStotal_time=cputime-ttime;


% fprintf('Space dim %d  Solution rank %d residual %d  time %d \n',size(VV,2),is,nrmres,RKStotal_time);


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)
 
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfuni(x,eH,s)

for j=1:length(x)
r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
    sval=linspace(eHpoints(j),eHpoints(j+1),40);
   %sval=linspace(eHpoints(j),eHpoints(j+1),100);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RA=compute_RA(A,eH,s);

ne=length(eH);
s=[1e20,s];
I=speye(size(A));
RA=I;
for k=1:ne,
   RA = (A-eH(k)*I)/(A-s(k)*I)*RA;
end
return

