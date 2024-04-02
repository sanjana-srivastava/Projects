% g(1), g(2)… -> inequality constraints
% h(1), h(2), …-> equality constraints
%
function [h,g] = constr(x)
h(1) = 0;
load = 6000;
length = 14;
modulusE = 30e6;
modulusG = 12e6;
tmax = 13600;
sigmamax = 30000;
delmax = 0.25;
tdash = load/(sqrt(2)*x(1)*x(2));
R = sqrt(x(2)*x(2)/4 + ((x(1)+x(3))/2)^2);
M = load*(length + x(2)/2);
J = 2* ((x(1)*x(2)/sqrt(2)) * (x(2)^2/12 +((x(1)+x(3))/2)^2));
tdashdash = M*R/J;
tx = sqrt(tdash^2 + 2*tdash*tdashdash*x(2)/(2*R) + tdashdash^2);
sigmax = 6*load*length/(x(4)*x(3)^2);
delx = 4*load*length^3/(modulusE*x(4)*x(3)^3);
pcx = (4.013*sqrt(modulusE*modulusG*x(3)^2*x(4)^6/36)/(length^2)) *...
(1- (x(3)/(2*length))*sqrt(modulusE/(4*modulusG)));
g(1) = tx/tmax -1;
g(2) = sigmax/sigmamax -1;
g(3) = x(1) - x(4);
g(4) = (.10471*x(1)*x(1) + 0.04811*x(3)*x(4)*(14+x(2)))/5 -1;
g(5) = 0.125 - x(1);
g(6) = delx/delmax -1;
g(7) = load/pcx -1;
g(8) = x(1)/2 -1;
g(9) = x(4)/2 -1;
g(10) = -x(1) + 0.1;
g(11) = -x(4) + 0.1;
g(12) = x(2)/10 -1;
g(13) = x(3)/10 -1;
g(14) = -x(2) + 0.1;
g(15) = -x(3) + 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%