%% Optimal Spacecraft Trajectories Project
% Time optimal Earth-Saturn transfer trajectory
close all; clear; clc;

%% Initialization

opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-14); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'TolX',1e-14,'Algorithm', 'levenberg-marquardt'); % fsolve

% Constants
AU = 149597870691;
TU = 58.13*24*60*60;
mu = 1; % gravity parameter
r_e = 1;
r_f = 9.5;
b_2 = 1;
P_re = 4.563e-6;
m = 300;
i = 0;

figure
polarplot(0,'o')
hold on
polarplot(0:0.01:2*pi,1*ones(1,length(0:0.01:2*pi)))
polarplot(0:0.01:2*pi,r_f*ones(1,length(0:0.01:2*pi)))
title('Earth-Saturn Transfer | Time Optimal Trajectories')

Amin = 100000;
A_increment = 20000;
Amax = 400000;

for A = Amin:A_increment:Amax

    i = i+1;
       
    a_c = ((2*P_re*A)/m)*((TU^2)/(AU));
    beta = a_c*r_e*b_2;
    r_f = 9.5;

    % Initial and final values
    x0 = [1 0 0 1.25];
    xf = [r_f 0 sqrt(mu/(r_f))];

    % time
    t0 = 0;

    % Minimize the time
    lam0_guess_tf = [-1 0 -1 -1 10];

    [lam0,fval] = fsolve(@cost_minT, lam0_guess_tf, options, t0, x0, xf, mu, beta, opts_ode);
    tf   = lam0(5);
    lam0 = lam0(1:4);
    [t,xl] = ode45(@(t,xl) dyn_lam(t,xl,mu,beta), [t0 tf], [x0 lam0], opts_ode);
    H = hamil(t, xl, mu, beta);
    H_minT = H';
    H_minT = H_minT+1;
    Tf(i) = t(end);
    
    lamr0 = xl(end,5);
    lamt0 = xl(end,6);
    lamu0 = xl(end,7);
    lamv0 = xl(end,8);

    polarplot(xl(:,2),xl(:,1),'LineWidth',1)

end

figure
Tf = Tf*TU; %  seconds
Tf = Tf /60; % minutes
Tf = Tf / 60; %hours
Tf = Tf / 24; %days
Tf = Tf / 365; % years
A = Amin:A_increment:Amax;
plot(A,Tf,'LineWidth',2,'Marker','*')
xlabel('Area (m^2)')
ylabel('Time of flight (years)')
title('Area effects of time of flight')
%% Hamiltonian
function H = hamil(t,xl,mu,beta)
    
    alp = atan(((-3*xl(:,7)) - sqrt(9*((xl(:,7)).^2) + 8*((xl(:,8)).^2)))./(4*xl(:,8)));
    xdot = [(xl(:,3)), (xl(:,4)./xl(:,1)), ((((xl(:,4)).^2)./(xl(:,1))) - (mu./((xl(:,1)).^2)) + (beta*((cos(alp)).^3))./(xl(:,1).^2)), (- (((xl(:,3)).*(xl(:,4)))./(xl(:,1))) + (beta*sin(alp).*((cos(alp)).^2))./(xl(:,1).^2))];
    lam = [xl(:,5) xl(:,6) xl(:,7) xl(:,8)];
    
    H = lam*xdot';
end

%% Function to minimize

function e = cost_minT(lam0_guess_tf,t0,x0,xf,mu,beta,opts_ode)
    lam0 = lam0_guess_tf(1:4);
    tf   = lam0_guess_tf(5);
    
    [t,xl] = ode45(@(t,xl) dyn_lam(t,xl,mu,beta), [t0 tf], [x0 lam0], opts_ode);
    H = hamil(t,xl,mu,beta);
    
    H = H';
    e = [xl(end,[1,3:4]) H(end)+1] - [xf 0];
end

%% Equation of differential of state and costate variables wrt time

function dxldt = dyn_lam(t,xl,mu,beta)
    alp = atan(((-3*xl(7)) - sqrt(9*((xl(7))^2) + 8*((xl(8))^2)))/(4*xl(8)));
   
    dxldt = zeros(8,1);
    % differentiation of the radial distance wrt time
    dxldt(1) = xl(3);
    % differentiation of the angular position wrt time
    dxldt(2) = xl(4)/xl(1);
    % differentiation of the radial velocity wrt time
    dxldt(3) = (((xl(4))^2)/(xl(1))) - (mu/((xl(1))^2)) + (beta*((cos(alp))^3))/(xl(1)^2);
    % differentiation of the tangential velocity wrt time
    dxldt(4) = - (((xl(3))*(xl(4)))/(xl(1))) + (beta*sin(alp)*((cos(alp))^2))/(xl(1)^2);
    % differentiation of lambda corresponding to radial distance wrt time
    dxldt(5) = (xl(7)/((xl(1))^3))*(2*beta*((cos(alp))^3) + (xl(4)^2)*xl(1) - 2*mu) + (xl(8)/((xl(1))^3))*(2*beta*sin(alp)*((cos(alp))^2) - xl(3)*xl(4)*xl(1));
    % differentiation of lambda corresponding to angular position wrt time
    dxldt(6) = 0;
    % differentiation of lambda corresponding to radial velocity wrt time
    dxldt(7) = - (xl(5) - ((xl(4)*xl(8))/(xl(1))));
    % differentiation of lambda corresponding to tangential velocity wrt time
    dxldt(8) =  ((- 2*xl(4)*xl(7) + xl(3)*xl(8))/(xl(1)));
end