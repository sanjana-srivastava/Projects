function [] = Trilobe_Design_Code()

clear                                                                      % clear all variable/information in the workspace - use CAUTION
clear global                                                               % again use caution - clears global information
clc                                                                        % position the cursor at the top of the screen
close                                                                      % closes the figure window
format compact                                                             % avoid skipping a line when writing to the command window
warning off                                                                % don't report any warnings like divide by zero etc.


timestamp = clock;
disp('Trilobed Airship Sizing Program');
disp(['Date: ',date,' Time: ',num2str(timestamp(4)),':', num2str(timestamp(5))]);
disp(' ');

%% start time for performance measure

starttime = cputime;

global Lmin Lmax Amin Amax hmin hmax IDAY XLONG XLAT YAW M_pay ...
       Objectivefunction Shapeselection Cdlift_yes noe Rho_fc
global Rho_batt ETA_SC P_pay k_pur Rho_en Solver N_lobes Filename 
global Destination storage

disp(' ');
 disp('Select shape of envelope:');
 disp(' ');
 disp(' 1: Ellipsoid');
 disp(' 2: NPL');
 disp(' 3: General');
 disp(' ');
 Shapeselection =input('Enter your selection: ');
 disp(' ');

Cdlift_yes = 1;                                                            % drag due to lift to be added (1) or not (0)
Lmin    = 20;                                                              % Minimum length (m)
Lmax    = 500;                                                             % Maximum length (m) [User defined - Lmin]
Amin    = 0;                                                               % Minimum angle of array from center (deg)
Amax    = 90;                                                              % Maximum angle of array from center (deg)
N_lobes = 3;                                                               % Number of lobes of envelope
noe     = 101;                                                             % Number of elements for solar array

disp(' ');
 disp('Select the city of operation:');
 disp(' ');
 disp(' 1: Mumbai');
 disp(' 2: Delhi');
 disp(' 3: Kolkata');
 disp(' 4: Chennai');
 disp(' ');
 city = input('The city to be operated: ');
 disp(' ');

switch city
    case 1
        XLONG = 72.88;                                                     % Longitude of operating location (deg)
        XLAT  = 19.07;                                                     % Latitude of operating location (deg)
        disp(('The selected operating location is Mumbai, India'));        % Display the detail of operating location

    case 2
        XLONG = 77.23;                                                     % Longitude of operating location (deg)
        XLAT  = 28.65;                                                     % Latitude of operating location (deg)
        disp(('The selected operating location is Delhi, India'));         % Display the detail of operating location
        
    case 3
        XLONG = 88.36;                                                     % Longitude of operating location (deg)
        XLAT  = 22.57;                                                     % Latitude of operating location (deg)
        disp(('The selected operating location is Kolkata, India'));       % Display the detail of operating location
        
    case 4
        XLONG = 80.24;                                                     % Longitude of operating location (deg)
        XLAT  = 13.07;                                                     % Latitude of operating location (deg)
        disp(('The selected operating location is Chennai, India'));       % Display the detail of operating location
end

disp(' ');
 disp('Select the energy storage device: ');
 disp(' ');
 disp(' 1: Rechargeable Battery');
 disp(' 2: Fuel cell');
 disp(' ');
 storage = input('Enter your selection: ');
 disp(' ');

switch storage
    case 1
        disp(('The selected energy storage device is Rechargeable Battery'));
    case 2
        disp(('The selected energy storage device is Fuel cell'));
end

% User-defined inputs
        
YAW      = 90;                                                             % Airship Yaw angle (deg)
hmin     = 0;                                                              % Minimum Altitude (km)
hmax     = 1;                                                              % Maximum Altiude  (km)
IDAY     = 209;                                                            % Day of operation
ETA_SC   = 0.12;                                                           % Solar cell eff
P_pay    = 1000;                                                           % Power to Payload(W)
Rho_batt = 400;                                                            % Energy density of battery(Wh/kg)
Rho_fc   = 1000;                                                           % Energy density of fuel cell (Wh/kg)
k_pur    = 0.97;                                                           % Percentage purity of helium (97%)
Rho_en   = 0.2;                                                            % Specific mass of envelope material - vectran(kg/m^2)

Dt    = datetime('1-Jan-2021')+IDAY-1;
disp(' ');
disp(['The selected day of operation is: ',datestr(Dt)]);
disp(' ');
M_pay = input('Desired Payload weight(kg): ');                             % Payload mass (kg)

disp(' ');
 disp('Select the objective function to be minimized: ');
 disp(' ');
 disp(' 1: Volume of the envelope');
 disp(' 2: Area of the array');
 disp(' 3: Mass of the airship');
 disp(' 4: Surface area of the envelope');
 disp(' 5: Drag coefficient of the envelope');
 disp(' ');
 Objectivefunction = input('Enter your selection: ');
 disp(' ');

% OptimizeTrilobe sends the requirement values to the optimization
% routine to minimize the desired objective.

disp(' ');
 disp('Select the optimizer/solver for optimization: ');
 disp(' ');
 disp(' 1: Particle Swarm Optimization (PSO)');
 disp(' 2: Genetic Algorithm (GA)');
 disp(' 3: Simulated Annealing (SIMANN)');
 disp(' 4: Grey Wolf Optimizer (GWO)');
 disp(' 5: Gaussian Quantum-based PSO (GQPSO)');
 disp(' ');
 
 Solver = input('Enter your selection: ');
 disp(' ');
 
 Filename = strcat('Airship');
 
 Destination = input('Name the Folder to store data: ','s');               % Folder name
 
%  if isempty(Desti)
%      Desti = Final_Output;
%  end

optimizeTrilobe()

%% save the optimized design variables in a text file

DesignVariable_Final()

%% postprocess displays all pertinent information about optimized design

postprocess()

%% print total time

totaltime = cputime - starttime;
fprintf('\n\nTotal cpu time (s)= %7.4f \n\n',totaltime)

end