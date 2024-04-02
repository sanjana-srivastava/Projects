% y -> objective function
% penalty -> penalty term
%
function y = func1(x)
y = 1.10471*x(1)*x(1)*x(2) + 0.04811*x(3)*x(4)*(14+x(2));
penalty = 0.0;
[h,g] = constr(x);
for i = 1:length(h)
if h(i)~=0
penalty = penalty + h(i)^2;
end
end
for i = 1:length(g)
if g(i)>0
penalty = penalty + g(i)^2;
end
end
scale_factor = 10000000;
y = y+penalty*scale_factor;