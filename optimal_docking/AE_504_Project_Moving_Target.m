%%% Optimal Docking Control AE 504 Project
%%% Madson Pelster Srivastava
%%% dmadson2
%%% AE504
%%% Dr. Negar Mehr

% The problem involves using thruster control to move the docking
% port and the space craft to the same location. The axis of
% the sc and dp must also align. For added complexity the target is
% non-stationary and moving randomly


% The spacecraft is subject to various constraints commented below and
% commented in the final report



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
rotation_lim = deg2rad(15) * 1/dt; % limits roation rate
docking_speed_lim = .5; % meters / second




%% Apply LQR Control function in various scenarios

% IN THIS CODE I WILL RE-RUN THE INITIAL SIMPLE TEST FUNCTION FOR DEBUG
% PURPOSES SINCE IT IS LOW COST

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

% CREATE A RANDOMLY MOVING TARGET IN 3 DIMESIONS


% CREATE A RANDOMLY MOVING AND ROTATING DOCKING PORT
% WE THEN SELECT THE MINIMUM TIME COST
random_rate_lim =  deg2rad(15) * 1/dt * .5; % target rotation rate lim


th = linspace(0,2*pi) ; 
R = 20 ;  % Radius 
C = [0 0 0] ; % center 
objx = C(1) + R*cos(th) ; 
objy = C(2) + R*sin(th) ;


objz = C(3) + 5*sin(th).* sin(th) ;

% plot(x,y)

% rotate the object randomly
obj_pos = [objx; objy; objz]; % position of target moving in a circle
obj_angle = zeros(3,100);
random1 = rand;
random2 = rand;
random3 = rand;

% get the random rotation rates to plug in later so we can match rotation
% rate1 = random1*random_rate_lim;
% rate2 = random2*random_rate_lim;
% rate3 = random3*random_rate_lim;

% for i = 1 : 99
% 
% obj_angle(1,i+1) = obj_angle(1,i) + random1*random_rate_lim;
% obj_angle(2,i+1) = obj_angle(2,i) + random2*random_rate_lim;
% obj_angle(3,i+1) = obj_angle(3,i) + random3*random_rate_lim;
% 
% 
% % make sure there are no values greater than 2*pi
% 
%     if obj_angle(1,i+1) >= 2*pi
%     
%     obj_angle(1,i+1) = obj_angle(1,i+1) / (2*pi);
%     
%     end
%     
%     if obj_angle(2,i+1) >= 2*pi
%     
%     obj_angle(2,i+1) = obj_angle(2,i+1) / (2*pi);
%         
%     end
%     
%     if obj_angle(3,i+1) >= 2*pi
%     
%     obj_angle(3,i+1) = obj_angle(3,i+1) / (2*pi);
%         
%     end
%     
% end


% NOW CALCULATE EACH INTERCEPT POSSIBILITY AND COMPARE THEM
for i = 1:size(obj_pos,2)

% varaible for total possible distance when
% genreating random scenarios 
dist = 200;
angle = 90;

% actual stating location (SUGGEST MAKING THE STARTING LOCATION
% RANDOM AND FINSIHING AT 0. MAY MAKE DEBUGGING EASIER)
PF = [obj_pos(1,i),obj_pos(2,i),obj_pos(3,i), ... 
    obj_angle(1,i),obj_angle(2,i),obj_angle(3,i)];

% location of docking port (finishing lcoation) {i just made this random 
% for initial scenario testing and debugging} we can
% create specific / complex scenarios later
PS = [dist,dist,dist, angle, ...
    angle,angle];


% starting of spacecraft position for scenario
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

intercept_point = [objx(min_point);objy(min_point);objz(min_point)];
intercept_angle = [obj_angle(1,min_point),obj_angle(2,min_point) ...
    obj_angle(3,min_point)];

fprintf(['\n The optimal time to intercept for the circular moving' ...
    'docking port is %f \n'], min_val)
fprintf(['\n The optimal intercept point is' ...
    ' %f \n'], intercept_point)


% note this only 2d so we will enhance functionality and make 3d in futue
% we will also rotate the docking objective randomly and recalculate

% make a nice plot of intercept location


figure(get(gcf,'Number') + 1)
hold on
plot3(objx,objy,objz)
plot3(objx(min_point),objy(min_point),objz(min_point),'o')
plot3(PS(1),PS(2),PS(3),'o')
title('Moving docking port',...
    'fontsize',14)
xlabel('X (meters)','fontsize',14)
ylabel('Y (meters)','fontsize',14)
zlabel('Z (meters)','fontsize',14)