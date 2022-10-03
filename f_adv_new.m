 function [y,TT] = f_adv_new(PDE,x)
    
 
    
    y = [ (x>=0).*x] .* (PDE.Tpos*x) + ...
        [(x<0).*x] .* (PDE.Tneg*x);
    %y = y;
if nargout ==2
TT =    [ (x>=0).*x] .* (PDE.Tpos) + ...
        [(x<0).*x] .* (PDE.Tneg);  
    %TT = TT;
else
    TT = [];
end
%keyboard 
 end  