%%% Optimal Docking Control AE 504 Project
%%% Madson Pelster Srivastava
%%% dmadson2
%%% AE504
%%% Dr. Negar Mehr

% Start with code that solves the classical LQR problem
% of optimal thrust control for aligning a space craft with a 
% docking port

% The problem involves using thruster control to move the docking
% port and the space craft to the same location. The axis of
% the sc and dp must also align. 

% {Perhaps an improvement on the complexity of the project can be that the
% docking port is spinning and the sc must match the spin to dock.

% TO DO:

% NEED TO PROGRAM A WAY TO GENERATE RANDOM SCENARIOS. GOES WITH INCRESAING
% COMPLEXITY. SHOULD BE EZ-PZ...

% NEED TO FIND SOMETHNG THAT MAKES OUR PROBLEM A BIT MORE COMPLEX
% 
% SUGGESTION IS TO CALCULATE LQR A BUNCH OF TIMES ON A MOVING OBJECT
% AND SELECT THE RESULT THAT USES THE LEAST AMOUT OF TIME



%%
  clear all
% close all

% global run_num fileID
% 
% run_num = 1;
% 
% name = 'Run 2.txt';
% 
% fileID = fopen(name,run_num);
%% Declare Constants

initial_tol = 5; % initial tolerance before aligning with docking port

Q = eye(6); % set Q as identity matrix. 
R = eye(6); % start with identity matrix and we have 6 control inputs

dt = .1; % time interval desired to be used during LQR processing
N = 50; % Number of iterations to run during LQR

docking_tol = .05; % tol for docking to know when final stage is complete
align_tol = .008; % tol for alignment stage (1/2 degree in rads)
rotation_lim = deg2rad(15) * 1/dt; %limits roation rate
docking_speed_lim = .5;




%% Apply LQR Control function in various scenarios

% here we establish the scenario then call
% each of the three different stage functions
% that in turn call LQR_Thruster to perform LQR
% optimization

% varaible for total possible distance when
% genreating random scenarios (WILL NEED TO BE
% MORE COMPLEX IN FUTURE. JUST TEMP MEASURE)
dist = 200;
angle = 90;

% actual stating location (SUGGEST MAKING THE STARTING LOCATION
% RANDOM AND FINSIHING AT 0. MAY MAKE DEBUGGING EASIER)
PF = zeros(1,6);

% location of docking port (finishing lcoation) {i just made this random 
% for initial scenario testing and debugging} we can
% create specific / complex scenarios later
PS = [dist*rand,dist*rand,dist*rand, angle*rand, ...
    angle*rand,angle*rand];


% starting of spacecraft position for scenario
% NEED TO ADD SCENARIO GENERATOR HERE
v(1)= 1;
x_dot(1) = 1;
y_dot(1) = 1;
z_dot(1) = 1;
alp = PS(4);
bet = PS(5);
gamm = PS(6);
% x(:,:,1) = PS;

% Write initial conditions to text document
% fprintf(fileID,'INITIAL STARTING CONDITIONS \n');
% fprintf(fileID,'%f ',PS);
% fprintf(fileID,'\n ',PS);
% fprintf(fileID,'Docking Port Location \n');
% fprintf(fileID,'%f ',PF);
% fprintf(fileID,'\n ',PF);

% Call initial stage for docking
[initial(1,1:3,1),tci] = LQR_Docking_Initial(dt,N,PS,PF,x_dot,...
    y_dot,z_dot,initial_tol);

time_initial = tci * dt; % gets total time for initial portion
% Call alignment stage for docking
[align(1,1:3,1),tca] = LQR_Docking_Align(dt,N,PS,PF,...
    alp,bet,gamm,align_tol,rotation_lim);

time_align = tca * dt;

% combine elements from above
x(1,1:3,1) = initial;
x(1,4:6,1) = align;

% call final stage with 6x6 LQR for docking
[x(1,:,end),tcf]= LQR_Docking_Final(dt,N,x(1,:,1),PF,...
    x_dot,y_dot,z_dot,x(1,4,1),...
    x(1,5,1),x(1,6,1),v,docking_tol,docking_speed_lim);

time_final = tcf * dt;

total_time = time_initial + time_align +time_final;

% Output the final resulting data:

pos_string = ['x_pos ','y_pos ','z_pos ','alpha ','beta ','gamma '];


fprintf('Final Docking Position Reached \n \n')
fprintf(['LQR Optimixation Control resulted in a final docking position ' ...
    'of \n \n:'])
fprintf('%g  ',x)



% fclose(fileID);




% To do: 

% add random pertubations to increase complexity

% create some basic scenarios and calculate results

% find a way to make it a bit more dynamical / interesting
% i.e how is this different than a hw problem


% random rotating target moving in a circle

th = linspace(0,2*pi) ; 
R = 20 ;  % Radius 
C = [0 0] ; % center 
objx = C(1)+R*cos(th) ; 
objy = C(2)+R*sin(th) ;
objz = zeros(100);

% plot(x,y)

obj_pos = [objx; objy; objz]; % position of target moving in a circle
obj_angle =[0,0,0];

% now run calculation for each 100 possible intercepts
for i = 1:size(obj_pos,2)

% varaible for total possible distance when
% genreating random scenarios (WILL NEED TO BE
% MORE COMPLEX IN FUTURE. JUST TEMP MEASURE)
dist = 200;
angle = 90;

% actual stating location (SUGGEST MAKING THE STARTING LOCATION
% RANDOM AND FINSIHING AT 0. MAY MAKE DEBUGGING EASIER)
PF = [obj_pos(1,i),obj_pos(2,i),obj_pos(3,i),obj_angle(1),obj_angle(2),obj_angle(3)];

% location of docking port (finishing lcoation) {i just made this random 
% for initial scenario testing and debugging} we can
% create specific / complex scenarios later
PS = [dist,dist,dist, angle, ...
    angle,angle];


% starting of spacecraft position for scenario
% NEED TO ADD SCENARIO GENERATOR HERE
v(1)= 1;
x_dot(1) = 1;
y_dot(1) = 1;
z_dot(1) = 1;
alp = PS(4);
bet = PS(5);
gamm = PS(6);
% x(:,:,1) = PS;

% Write initial conditions to text document
% fprintf(fileID,'INITIAL STARTING CONDITIONS \n');
% fprintf(fileID,'%f ',PS);
% fprintf(fileID,'\n ',PS);
% fprintf(fileID,'Docking Port Location \n');
% fprintf(fileID,'%f ',PF);
% fprintf(fileID,'\n ',PF);

% Call initial stage for docking
[initial(1,1:3,1),tci] = LQR_Docking_Initial(dt,N,PS,PF,x_dot,...
    y_dot,z_dot,initial_tol);

time_initial = tci * dt; % gets total time for initial portion
% Call alignment stage for docking
[align(1,1:3,1),tca] = LQR_Docking_Align(dt,N,PS,PF,...
    alp,bet,gamm,align_tol,rotation_lim);

time_align = tca * dt;

% combine elements from above
x(1,1:3,1) = initial;
x(1,4:6,1) = align;

% call final stage with 6x6 LQR for docking
[x(1,:,end),tcf]= LQR_Docking_Final(dt,N,x(1,:,1),PF,...
    x_dot,y_dot,z_dot,x(1,4,1),...
    x(1,5,1),x(1,6,1),v,docking_tol,docking_speed_lim);

time_final = tcf * dt;

total_time_moving(i) = time_initial + time_align +time_final;


end

[min_val, min_point] = min(total_time_moving);

intercept_point = [objx(min_point);objy(min_point)];

fprintf(['\n The optimal time to intercept for the circular moving' ...
    'docking port is %f \n'], min_val)
fprintf(['\n The optimal intercept point is' ...
    ' %f \n'], intercept_point)


% note this only 2d so we will enhance functionality and make 3d in futue
% we will also rotate the docking objective randomly and recalculate

% make a nice plot of intercept location

figure(2)
hold on
title('LQR Optimal Intercept for Circular Moving Docking Port',...
    'fontsize',14)
plot(objx,objy)
plot(obj_pos(1,27), obj_pos(2,27),'o')
plot(PS(1),PS(2),'o')
xlabel('X (meters)','fontsize',14)
ylabel('Y (meters)','fontsize',14)
legend({'Docking Port Path','Intercept Position',...
    'Start Location'},'Location','southeast')










