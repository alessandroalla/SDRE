function  [Tpos,Tneg] = get_adv(n)
%global Txup Tyup Txdown Tydown
dx = 1/(n+1);
dy = dx;
uno = ones(n,1);
Tpos = eye(n) - diag(ones(n-1,1),-1); % se vel pos
Tpos = sparse(Tpos);
Tneg = diag(ones(n-1,1),1) - eye(n);
Tneg = sparse(Tneg);
%Tpos(1,n) = -Tpos(1,1);
%Tneg(n,1) = -Tneg(n,n);
Tpos = Tpos/dx;
Tneg = Tneg/dx;
Tpos = kron(Tpos,speye(n)) + kron(speye(n),Tpos);
Tneg = kron(Tneg,speye(n)) + kron(speye(n),Tneg);

% Txup = sparse(n^2,n^2);
%   Tyup = Txup;
%   Txdown=Txup;
%   Tydown=Txup;
%   for ix = 1:n
%     for iy = 1:n
%        % Entry for x-direction:
%  %      keyboard
%         idx = (iy-1)*n + ix;
%          Txup(idx,idx) = 1/dx;
%                     Txdown(idx,idx) = -1/dx;
%                     if ix > 1
%                         Txup(idx,idx-1) = -1/dx;
%                     end
%                     if ix < n
%                         Txdown(idx,idx+1) = 1/dx;
%                     end
%                     
%                     % Entry for y-direction:
%                     idx = (iy-1)*n+ix;
%                     Tyup(idx,idx) =  1/dx;
%                     Tydown(idx,idx) = -1/dx;
%                     if iy > 1
%                         Tyup(idx,idx-n) = -1/dx;
%                     end
%                     if iy < n
%                         Tydown(idx,idx+n) = 1/dx;
%                     end
%                 end
%             end
%             Txup = Txup;
%             Tyup = Tyup;
%             Txdown = Txdown;
%             Tydown = Tydown;
% %            Txdown(1,n^2) = - Txdown(1,1);
% %            Txup(n^2,1) = -Txup(n^2,n^2);
% %            Txup(1,n^2) = - Txup(1,1);
% %            Txdown(n^2,1) = -Txdown(n^2,n^2);
            

end

