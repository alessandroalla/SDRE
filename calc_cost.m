function [cost] = calc_cost(PDE,y)

value(1) = ((PDE.C_tmp*y(:,1))'*(PDE.C_tmp*y(:,1)));
for i =2:size(y,2)
    value(i) =  ((PDE.C_tmp*y(:,i))'*(PDE.C_tmp*y(:,i)));
end

cost = sum(value);

end

