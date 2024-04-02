
%%% Optimal Docking Control AE 504 Project
%%% Madson Pelster Srivastava
%%% AE504
%%% Dr. Negar Mehr

% LQR ALIGNMENT CONTROLLER FOR DOCKING
% ALIGN SPACECRAFT WITH DOCKING PORT 
% WHEN WITHIN ACCEPTABLE DISTANCE
% USING REACTION THRUSTERS OR INERTIAL CONTROL SYSTEM


function [out,tc] = LQR_Docking_Align(dt,N,PS,PF,...
    alp,bet,gamm,align_tol,rotation_lim)

% constants based on 3d thrusting ablity (seek to align angle only
% not position)

tc = 0;

 Q = eye(3);
 R = eye(3);

 x = nan(1,3,100);
 x(:,:,1) = PS(4:6);
 PS(1:3) = PS(4:6);


i = 1; % counter to determine what time step we are on

% THE BELOW WHILE LOOPS RUNS UNTIL THE SPACECRAFT IS IN THE
% CORRECT POSTION FOR DOCKING. WE CAN CHANGE THE TOLERANCE LEVEL LATER

check = norm(PS(1:3)) - norm(PF(4:6));
while abs( PF(4) - x(1,1,i) ) > align_tol || abs( PF(5) - x(1,2,i)) > align_tol || abs( PF(6) - x(1,3,i) ) > align_tol

    i= i+1;
   
% A matrix and B matrix
% here we need to calculate our A and B matrices
% these are essentially the equations of motion of
% the space craft.

A = [1 0 0;...
    0 1 0;...
    0 0 1];

B = [dt 0 0;...
    0 dt 0;...
    0 0 dt];



% Run and check the LQR_Thruster function to see what the
% next control input should be
% 
% fprintf('making LQR decision for alignment stage %g',(i-1))

move = LQR_Thruster(A,B,Q,R,dt,PS(1:3),PF(4:6),N);

us = move;

% Establish a rotation limit

    limit = us(1:3);
    limit(limit > rotation_lim) = rotation_lim;
    limit(limit < -rotation_lim) = -rotation_lim;
    us(1:3) = limit;

% move the spacecraft based on resulted u_star
% we will apply the state space model

% x_t = A(t-1)*x(t-1)+B(t-1)u(t-1)

x(:,:,i) = A*x(:,:,i-1)' + B*us;
% v(i) = us(1); %??? not sure if this is right

% update positon based on caclulated movements
PS = x(:,:,i);

% set PS (can't do this right now until I make more progress on what
% u* or us will look like when it comes out and if the math we have even
% makes sense)

check = norm(PS) - norm(PF);

tc = tc + 1;

end

out = x(:,:,i);






end










