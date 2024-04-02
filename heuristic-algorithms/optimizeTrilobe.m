function []= optimizeTrilobe()

global u hmin hmax fbest Shapeselection Solver


% Design Variables:-

% m   : Ratio of distance of max.section from the nose to length
% R0  : Radius of curvature at nose of hull
% R1  : Radius of curvature at tail of hull
% Cp  : Prismatic coefficient
% L/D : Fineness ratio
% f   : Relative distance (Lateral)
% g   : Relative distance (Vertical)

switch Shapeselection
    case 1
        Mmin   = 0.500; 
        Mmax   = 0.500;
        r0min  = 0.500;
        r0max  = 0.500;
        r1min  = 0.500;
        r1max  = 0.500;
        cpmin  = 0.666;
        cpmax  = 0.666;
        l2dmin = 4.999;
        l2dmax = 4.999;
    case 2
        Mmin   = 0.431939; 
        Mmax   = 0.431939;
        r0min  = 0.5886;
        r0max  = 0.5886;
        r1min  = 0.42485;
        r1max  = 0.42485;
        cpmin  = 0.66675;
        cpmax  = 0.66675;
        l2dmin = 4;
        l2dmax = 4;
    case 3
        Mmin   = 0.3; 
        Mmax   = 0.6;
        r0min  = 0.5;
        r0max  = 1.0;
        r1min  = 0;
        r1max  = 0.5;
        cpmin  = 0.55;
        cpmax  = 0.70;
        l2dmin = 3;
        l2dmax = 6;
end

u6min  = 0;
u6max  = 1;
u7min  = 0;
u7max  = 1;
u8min  = 0;
u8max  = 1;
u10min = 0;
u10max = 1;
fmin   = 0.1;
fmax   = 0.3;
gmin   = 0.0;
gmax   = 0.2;
w1min  = 0.0;
w1max  = 0.9;
w2min  = 0.3;
w2max  = 0.3;
w3min  = 0.0;
w3max  = 0.9;
w4min  = 0.0;
w4max  = 0.9;

%% Design Parameters Bound

lb = [Mmin, r0min, r1min, cpmin, l2dmin, u6min, u7min, u8min, hmin, u10min, fmin, gmin, w1min, w2min, w3min, w4min];
ub = [Mmax, r0max, r1max, cpmax, l2dmax, u6max, u7max, u8max, hmax, u10max, fmax, gmax, w1max, w2max, w3max, w4max];

%% Number of design variables

nvars = length(lb);

switch Solver
    case 1
        % PSO
        rng(3, 'twister');
        hybridopts = optimoptions('patternsearch','Display','iter',...
                                  'InitialMeshSize',nvars*500,...
                                  'ScaleMesh',false,...
                                  'AccelerateMesh',true,...
                                  'MeshTolerance', 1e-7);
%                                 'PlotFcn',{@psplotbestf,@psplotfuncount});
% For fmincon

% hybridopts = optimset('Algorithm','interior-point','MaxFunEvals',5000,'TolCon',.01,'TolX',1e-9);
                      
%% Options for PSO

options = optimoptions('particleswarm','Display','iter',...
                       'SwarmSize',nvars*10,...
                       'MaxIterations',100,...
                       'PlotFcn',@pswplotbestf,...
                       'MinNeighborsFraction',1,...
                       'SelfAdjustmentWeight',0,...
                       'SocialAdjustmentWeight',0,...
                       'UseVectorized', false,...
                       'HybridFcn',{@patternsearch,hybridopts});


[xbest, fbest,~,~] = particleswarm(@Optim_Trilobe,nvars,lb,ub,options);

% [xbest,fbest,ExitFlag]=fmincon(@Optim_Trilobe,x_0,[],[],[],[],lb,ub,[],options);

    case 2
        % Options for GA
        hybridopts = optimoptions('patternsearch','Display','iter',...
                                  'InitialMeshSize',nvars*500,...
                                  'ScaleMesh',false,...
                                  'AccelerateMesh',true,...
                                  'MeshTolerance', 1e-7);
%                                 'PlotFcn',{@psplotbestf,@psplotfuncount});
        opts = gaoptimset(...
                         'Display','iter','PlotInterval',1,...
                         'PopulationSize', nvars*10, ...
                         'Generations',10,...
                         'EliteCount',[1], ...  
                         'CrossoverFraction', 0.9, ...
                         'MutationFcn', {@mutationuniform, 0.01},...
                         'TolFun', 1e-8, ...
                         'MigrationFraction',0.01, ...
                         'PlotFcns', @gaplotbestf, ...
                         'HybridFcn',{@patternsearch,hybridopts});
                     
rng(3, 'twister');
[xbest, fbest, ~] = ga(@Optim_Trilobe, nvars, [], [], [], [],lb,ub, [],[], opts);

    case 3
        [xbest] = simann('Optim_Trilobe', [0.6,0.6,0.6,0.6,6,0.3,0.7,0.3,0.5,0.4,0.4,0.4,0.6,0.6,0.6,0.9],lb,ub, 5000, .5, 5, 20, 100);
        feval('Optim_Trilobe',xbest)
        
    case 4
        SearchAgents_no = 30;                                              % Number of search agents
        fobj = @Optim_Trilobe;                                             % Name of the test function
        Max_iteration = 500;                                               % Maximum numbef of iterations
        dim = nvars;                                                       % Number of design variables
        [fbest,xbest,~] = GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    case 5
       [fbest, xbest] = GAQPSO(@Optim_Trilobe);
end

%%

x(1)  = xbest(1);
x(2)  = xbest(2);
x(3)  = xbest(3);
x(4)  = xbest(4);
x(5)  = xbest(5);
x(6)  = xbest(6);
x(7)  = xbest(7);
x(8)  = xbest(8);
x(9)  = xbest(9);
x(10) = xbest(10);
x(11) = xbest(11);
x(12) = xbest(12);
x(13) = xbest(13);
x(14) = xbest(14);
x(15) = xbest(15);
x(16) = xbest(16);

u = [x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9) x(10) x(11) x(12) x(13) x(14) x(15) x(16)];

end
