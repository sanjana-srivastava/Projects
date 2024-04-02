%%% Optimal Docking Control AE 504 Project
%%% Madson Pelster Srivastava
%%% AE504
%%% Dr. Negar Mehr

% LQR INITIAL CONTROLLER FOR DOCKING
% MOVES SPACECRAFT TO ACCURATE DOCKING DISTANCE
% USING 3-D RCS THRUSTERS


function [out,tc] = LQR_Docking_Initial(dt,N,PS,PF,x_dot,y_dot,z_dot,initial_tol)

% fprintf(fileID,'BEGENNING DOCKING INITIAL STAGE \n');


% constants based on 3d thrusting ablity (seek to move position only
% not align angle)

 tc = 0;

 Q = eye(3);
%  Q = [1 0 0; 0 10 0; 0 0 10];
 R = eye(3);

 x = nan(1,3,100);
 x(:,:,1) = PS(1:3);


i = 1; % counter to determine what time step we are on

% THE BELOW WHILE LOOPS RUNS UNTIL THE SPACECRAFT IS IN THE
% CORRECT POSTION FOR DOCKING. WE CAN CHANGE THE TOLERANCE LEVEL LATER

check = norm(PS(1:3)) - norm(PF(1:3));
while abs( PF(2) - x(1,2,i) ) > .25 || abs( PF(3) - x(1,3,i) ) > .25
%   x(1,:,i)
    i= i+1;
   

% fprintf(fileID,'BEGENNING DOCKING INITIAL STAGE \n');

% A matrix and B matrix
% here we need to calculate our A and B matrices
% these are essentially the equations of motion of
% the space craft.

A = [1 0 0;...
    0 1 0;...
    0 0 1];

% THE CODE BELOW STOPS THE MOTION ON THE X AXIS AND SHIFTS TO 
% ALIGNMENT VIA THE MAIN CODE WHEN THE SPACECRAFT IS WITHIN THE TOLERANCE
if abs(PF(1) - x(1,1,i-1) ) < (initial_tol + 40)

Q = [1 0 0; 0 10 0; 0 0 10];

end

if abs(PF(1) - x(1,1,i-1) ) < initial_tol

B = [0 0 0;...
    0 y_dot*dt 0;...
    0 0 z_dot*dt];



else

    B = [x_dot*dt 0 0;...
    0 y_dot*dt 0;...
    0 0 z_dot*dt];
    

end
    



% Run and check the LQR_Thruster function to see what the
% next control input should be

% fprintf(fileID,'Making LQR decision initial stage step %f \n',(i-1))


move = LQR_Thruster(A,B,Q,R,dt,PS(1:3),PF(1:3),N);


us = move;

% THE CODE BELOW INSTALLS A SPEED LIMIT WHEN THE SPACE CRAFT IS 
% WITHIN 20 METERS ON EITHER THE X Y OR Z AXIS



if check > 20
    
    while norm(us) / dt > 40
        us = us ./ 1.5;
    end
%     limit = us(1:3);
%     limit(limit > 5) = 5;
%     limit(limit < -5) = -5;
%     us(1:3) = limit;


else
    
    while norm(us) / dt > 5
        us = us ./ 1.5;
    end
    
%     limit = us(1:3);
%     limit(limit > 1) = 1;
%     limit(limit < -1) = -1;
%     us(1:3) = limit;
end


% move the spacecraft based on resulted u_star
% we will apply the state space model

% x_t = A(t-1)*x(t-1)+B(t-1)u(t-1)

x(:,:,i) = A*x(:,:,i-1)' + B*us;
v(i) = us(1); %??? not sure if this is right

% fprintf(fileID,'NEW POSITION AFTER STEP %f is: \n',(i-1))
% fprintf(fileID,'%f  ',x(:,:,i))
% fprintf(fileID,'\n ',x)

PS = x(:,:,i);

tc = tc + 1; %time count for export to main function

% set PS (can't do this right now until I make more progress on what
% u* or us will look like when it comes out and if the math we have even
% makes sense)

check = norm(PS) - norm(PF);


end

out = x(:,:,i);




for j = 1:i

    plot_x(j) = x(1,1,j);
    plot_y(j) = x(1,2,j);
    plot_z(j) = x(1,3,j);
    



end

    figure(get(gcf,'Number') + 1)
    plot3(plot_x,plot_y,plot_z)
    axis([0 200 0 200 0 200])
    title('LQR Optimal Intercept for Circular Moving Docking Port',...
    'fontsize',14)
    xlabel('X (meters)','fontsize',14)
    ylabel('Y (meters)','fontsize',14)
    zlabel('Z (meters)','fontsize',14)
%     






end












