function [jaco] = jacobiano(f,y,tmp,K,PDE,t,g,v,app)
%function [jaco] = jacobiano(f,y,tmp,K,PDE,t,g,v,app)

par = 1e-6;
jaco = f(y,tmp+par*v,K,PDE,t,g) - app;

jaco = jaco/(par);
end

