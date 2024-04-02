%% Optimal Spacecraft Trajectories HW 3
% Minimize control problem
close all; clear; clc;

%% Initialization
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-14); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-14,'Algorithm', 'levenberg-marquardt'); % fsolve

% Constants
mu = 1; % gravity parameter

% Initial and final values
x0 = [1 0 0 1.25];
xf = [9.5 0 sqrt(mu/(9.5))];

% time
t0 = 0;
tf_days = 1800;
tf = tf_days/58.13;

%% Minimize the control parameter

lam0_guess = [-0.01 -0.01 -0.01 -0.01];

[lam0,fval] = fsolve(@cost_minJ,lam0_guess,options,t0,tf,x0,xf);

[t,xmin] = ode45(@(t,xl) dyn_lam(t,xl), [t0 tf], [x0 lam0], opts_ode);

%% accumulative cost

L = 0.5*(xmin(:,7).^2+xmin(:,8).^2);
J = cumtrapz(L);

%% Plotting

% Trajectory plot
figure
polarplot(xmin(:,2),xmin(:,1),'LineWidth',2)
hold on
polarplot(0,'o')
polarplot(0:0.01:2*pi,1*ones(1,length(0:0.01:2*pi)))
polarplot(0:0.01:2*pi,1.524*ones(1,length(0:0.01:2*pi)))
title('Earth-Mars Transfer | Control Optimal Trajectory')
legend("Optimal Trajectory","Sun","Earth's Orbit","Mars's Orbit")

% States
figure
subplot(2,2,1)
plot(t,xmin(:,1),LineWidth=2)
xlabel('t')
ylabel('r')
title('Radial distance vs. Time')

subplot(2,2,2)
plot(t,xmin(:,2)*180/pi,LineWidth=2)
xlabel('t')
ylabel('\theta in ^o')
title('Angular position vs. Time')

subplot(2,2,3)
plot(t,xmin(:,3),LineWidth=2)
xlabel('t')
ylabel('u')
title('Radial velocity vs. Time')

subplot(2,2,4)
plot(t,xmin(:,4),LineWidth=2)
xlabel('t')
ylabel('v')
title('Tangential velocity vs. Time')

% Costates
figure
subplot(2,2,1)
plot(t,xmin(:,5),LineWidth=2)
xlabel('t')
ylabel('\lambda_r')
title('\lambda_r vs. Time')

subplot(2,2,2)
plot(t,xmin(:,6),LineWidth=2)
xlabel('t')
ylabel('\lambda_\theta')
title('\lambda_\theta vs. Time')

subplot(2,2,3)
plot(t,xmin(:,7),LineWidth=2)
xlabel('t')
ylabel('\lambda_u')
title('\lambda_u vs. Time')

subplot(2,2,4)
plot(t,xmin(:,8),LineWidth=2)
xlabel('t')
ylabel('\lambda_v')
title('\lambda_v vs. Time')

% Control
figure
subplot(2,1,1)
plot(t,-xmin(:,7),LineWidth=2)
xlabel('t')
ylabel('u_r')
title('u_r vs. Time')

subplot(2,1,2)
plot(t,-xmin(:,8),LineWidth=2)
xlabel('t')
ylabel('u_\theta')
title('u_\theta vs. Time')

% Accumulative cost
figure
plot(t,J,LineWidth=2)
xlabel('t')
ylabel('J')
title('Cost vs. Time')

%% Function to minimize
function e = cost_minJ(lam0,t0,tf,x0,xf)

    opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-14); % ode

    [t,xl] = ode45(@(t,xl) dyn_lam(t,xl), [t0 tf], [x0 lam0], opts_ode);

    e = xl(end,[1,3:4]) - xf;

end

%% Equation of differential of state and costate variables wrt time

function dxldt = dyn_lam(t,xl)

    mu = 1;
   
    dxldt = zeros(8,1);

    % differentiation of the radial distance wrt time
    dxldt(1) = xl(3);

    % differentiation of the angular position wrt time
    dxldt(2) = xl(4)/xl(1);

    % differentiation of the radial velocity wrt time
    dxldt(3) = (((xl(4))^2)/(xl(1))) - (mu/((xl(1))^2)) - xl(7);

    % differentiation of the tangential velocity wrt time
    dxldt(4) = - (((xl(3))*(xl(4)))/(xl(1))) - xl(8);

    % differentiation of lambda corresponding to radial distance wrt time
    dxldt(5) = (1/(xl(1)^2))*((xl(7)*xl(4)^2) - (2*mu*xl(7))/xl(1) - xl(8)*xl(3)*xl(4));

    % differentiation of lambda corresponding to angular position wrt time
    dxldt(6) = 0;

    % differentiation of lambda corresponding to radial velocity wrt time
    dxldt(7) = - (xl(5) - ((xl(4)*xl(8))/(xl(1))));

    % differentiation of lambda corresponding to tangential velocity wrt time
    dxldt(8) =  ((- 2*xl(4)*xl(7) + xl(3)*xl(8))/(xl(1)));

end