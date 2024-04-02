%% Optimal Spacecraft Trajectories Project
% Time optimal Earth-Saturn transfer trajectory with free \theta_f and ideal solar sail


%% Initialization
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-14); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'TolX',1e-14,'Algorithm', 'levenberg-marquardt'); % fsolve

% Constants
AU = 149597870691;
TU = 58.13*24*60*60;
mu = 1; % gravity parameter
r_e = 1;
b_2 = 1;
P_re = 4.563e-6;
m = 300;
A = 140000;

a_c = ((2*P_re*A)/m)*((TU^2)/(AU));

beta = a_c*r_e*b_2;

r_f = 9.5;

% Initial and final values
x0 = [1 0 0 1];
xf = [r_f 0 sqrt(mu/(r_f))];

% time
t0 = 0;

%% Minimize the time

lam0_guess_tf = [-1 0 -1 -1 10];

[lam0,fval] = fsolve(@cost_minT, lam0_guess_tf, options, t0, x0, xf, mu, beta, opts_ode);

tf   = lam0(5);
lam0 = lam0(1:4);

[t,xl] = ode45(@(t,xl) dyn_lam(t,xl,mu,beta), [t0 tf], [x0 lam0], opts_ode);

H = hamil(t, xl, mu, beta);
H_minT = H';

%% Plotting

% Trajectory plot
figure
polarplot(xl(:,2),xl(:,1),'LineWidth',2)
hold on
polarplot(0,'o')
polarplot(0:0.01:2*pi,1*ones(1,length(0:0.01:2*pi)))
polarplot(0:0.01:2*pi,r_f*ones(1,length(0:0.01:2*pi)))
title('Earth-Saturn Transfer | Time Optimal Trajectory')
legend("Optimal Trajectory","Sun","Earth's Orbit","Saturn's Orbit")

% States
figure
subplot(2,2,1)
plot(t,xl(:,1),LineWidth=2)
xlabel('t in TU')
ylabel('r in AU')
title('Radial distance vs. Time')

subplot(2,2,2)
plot(t,xl(:,2)*180/pi,LineWidth=2)
xlabel('t in TU')
ylabel('\theta in ^o')
title('Angular position vs. Time')

subplot(2,2,3)
plot(t,xl(:,3),LineWidth=2)
xlabel('t in TU')
ylabel('u in AU/TU')
title('Radial velocity vs. Time')

subplot(2,2,4)
plot(t,xl(:,4),LineWidth=2)
xlabel('t in TU')
ylabel('v in AU/TU')
title('Tangential velocity vs. Time')

% Costates
figure
subplot(2,2,1)
plot(t,xl(:,5),LineWidth=2)
xlabel('t in TU')
ylabel('\lambda_r')
title('\lambda_r vs. Time')

subplot(2,2,2)
plot(t,xl(:,6),LineWidth=2)
xlabel('t in TU')
ylabel('\lambda_\theta')
title('\lambda_\theta vs. Time')

subplot(2,2,3)
plot(t,xl(:,7),LineWidth=2)
xlabel('t in TU')
ylabel('\lambda_u')
title('\lambda_u vs. Time')

subplot(2,2,4)
plot(t,xl(:,8),LineWidth=2)
xlabel('t in TU')
ylabel('\lambda_v')
title('\lambda_v vs. Time')

% Control
alp = atan(((-3*xl(:,7)) - sqrt(9*((xl(:,7)).^2) + 8*((xl(:,8)).^2)))./(4*xl(:,8)));
figure
plot(t*(58.13/365.25),alp*180/pi,LineWidth=2)
xlabel('t in years')
ylabel('\alpha in ^o')
title('Pitch angle vs. Time')

% accelerations due to SRP
figure
subplot(2,1,1)
rasrp = (beta*((cos(alp)).^3))./(xl(:,1).^2);
plot(t*58.13,rasrp*(AU/(TU)^2),LineWidth=2)
xlabel('t in days')
ylabel('SRP a_r in m/s^2')
title('Radial acceleration due to SRP vs. Time')

subplot(2,1,2)
tasrp = (beta.*sin(alp).*((cos(alp)).^2))/(xl(1).^2);
plot(t*58.13,tasrp*(AU/(TU)^2),LineWidth=2)
xlabel('t in days')
ylabel('SRP a_t in m/s^2')
title('Tangential acceleration due to SRP vs. Time')

figure
resacc_srp = sqrt(rasrp.^2 + tasrp.^2);
plot3(t*(58.13/365.25),alp*180/pi,resacc_srp*(AU/TU^2),LineWidth=2)
hold on; grid on;
xlabel('t in years')
ylabel('\alpha in ^o')
zlabel('SRP a_R in m/(s^2)')
title('Resultant acceleration due to SRP vs. pitch angle')

% Hamiltonian
figure
subplot(2,1,1)
plot(t*58.13,H_minT,'r-','Linewidth',2)
ylabel('H')
title('Hamiltonian vs. time for time optimal case')

subplot(2,1,2)
semilogy(t*58.13,abs((H_minT-H_minT(1))/H_minT(1)),'r-','Linewidth',2)
hold on; grid on;
ylabel('|H-H(t_0)/H(t_0)|')
xlabel('t in days')

%% Hamiltonian

function H = hamil(t,xl,mu,beta)
    
    for i =1:length(xl)
        alp = atan(((-3*xl(i,7)) - sqrt(9*((xl(i,7))^2) + 8*((xl(i,8))^2)))/(4*xl(i,8)));

        xdot = [(xl(i,3)), (xl(i,4)/xl(i,1)), ((((xl(i,4))^2)/(xl(i,1))) - (mu/((xl(i,1))^2)) + (beta*((cos(alp))^3))/(xl(i,1)^2)), (- (((xl(i,3))*(xl(i,4)))/(xl(i,1))) + (beta*sin(alp)*((cos(alp))^2))/(xl(i,1)^2))];
        lam = [xl(i,5) xl(i,6) xl(i,7) xl(i,8)];
    
        H(i) = lam*xdot';
    end

end

%% Function to minimize

function e = cost_minT(lam0_guess_tf,t0,x0,xf,mu,beta,opts_ode)

    lam0 = lam0_guess_tf(1:4);
    tf   = lam0_guess_tf(5);
    
    [t,xl] = ode45(@(t,xl) dyn_lam(t,xl,mu,beta), [t0 tf], [x0 lam0], opts_ode);

    H = hamil(t,xl,mu,beta);
    
    H = H';

    e = xl(end,[1,3:4]) - xf;

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