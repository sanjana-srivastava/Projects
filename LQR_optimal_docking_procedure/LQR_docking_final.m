%%% Optimal Docking Control AE 504 Project
%%% Madson Pelster Srivastava
%%% AE504
%%% Dr. Negar Mehr

% LQR FINAL CONTROLLER FOR DOCKING
% MOVES SPACECRAFT TO DOCKING PORT FROM SHORT RANGE
% USING 3-D RCS THRUSTERS


function [out,tc] = LQR_Docking_Final(dt,N,PS,PF,...
    x_dot,y_dot,z_dot,alp,bet,gamm,v,docking_tol,docking_speed_lim)

tc = 0;

Q = eye(6);
R = eye(6);

x = nan(1,6,100);
x(1,:,1) = PS;



i=1;

check = norm(PS(1:3)) - norm(PF(1:3));
while check > docking_tol

    i= i+1;
   
% A matrix and B matrix
% here we need to calculate our A and B matrices
% these are essentially the equations of motion of
% the space craft.

A = [1 0 0 -v(i-1)*sin(alp)*dt 0 0;...
    0 1 0 0 -v(i-1)*sin(bet)*dt 0;...
    0 0 1 0 0 -v(i-1)*sin(gamm)*dt;...
    0 0 0 1 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1]; % the sin terms may be incorrect. Need to check
% need to make sure that this models how the space craft moves

B = [x_dot 0 0 0 0 0;...
    0 y_dot 0 0 0 0;...
    0 0 z_dot 0 0 0;...
    0 0 0 dt 0 0;...
    0 0 0 0 dt 0;...
    0 0 0 0 0 dt];
    



% Run and check the LQR_Thruster function to see what the
% next control input should be

% fprintf('making LQR decision step %g',(i-1));

move = LQR_Thruster(A,B,Q,R,dt,PS,PF,N);


us = move;

% slow the rate of motion as required near docking port
if check > 20
    limit = us(1:3);
    limit(limit > 5) = 5;
    limit(limit < -5) = -5;
    us(1:3) = limit;

else
    limit = us(1:3);
    limit(limit > docking_speed_lim) = docking_speed_lim;
    limit(limit < -docking_speed_lim) = -docking_speed_lim;
    us(1:3) = limit;


end


% move the spacecraft based on resulted u_star
% we will apply the state space model

% x_t = A(t-1)*x(t-1)+B(t-1)u(t-1)

x(:,:,i) = A*x(:,:,i-1)' + B*us;
v(i) = us(1);

alp = x(1,4,i);
bet = x(1,5,i);
gamm = x(1,6,i);


PS = x(:,:,i);

% set PS based on new caclulated position

check = norm(PS(1:3)) - norm(PF(1:3));

tc = tc + 1;

end

out = x(:,:,i);

end