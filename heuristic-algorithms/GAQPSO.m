%% management functions
function [BestFx, BestX]=GAQPSO(fun);
        % clear all variable/information in the workspace - use CAUTION   % again use caution - clears global information
         % position the cursor at the top of the screen
           % closes the figure window
format compact  % avoid skipping a line when writing to the command window
warning off     % don't report any warnings like divide by zero etc.

D = 4;                                                                     % Number of dimensions
nPop = 200;                                                                 % Population size
lb = [0.1,0.1,0.1,0.1];                                                    % Lower bounds
ub = [2,10,10,2];                                                          % Upper bounds 
itermax = 1500;                                                            % Maximum iterations
c1 = 2.05;                                                                 % Constant 1
c2 = 2.05;                                                                 % Constant 2
w1 = 1.0;
w2 = 1.0;

%% Constriction Coefficients
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2/(phi-2+sqrt(phi^2-4*phi));
w1 = chi;          % Inertia Weight
w2 = chi;
% wdamp=1;        % Inertia Weight Damping Ratio
c1 = chi*phi1;    % Personal Learning Coefficient
c2 = chi*phi2;    % Global Learning Coefficient


tic

%% Generating the initial population

pos = zeros(nPop,D);

for i = 1:nPop
    for j = 1:D
        pos(i,j) = lb(j)+rand.*(ub(j)-lb(j));
    end
end

%% Evaluate the function

fx    = fun(pos);

%% Initialize Pbest

pbest = pos;

%% Initialize Gbest

[fminvalue,ind] = min(fx);
gbest           = pos(ind,:);

%% Main Loop

iter = 1;
while iter < itermax
    Fgbest = fminvalue;
    beta = w1+(w2-w1)*(itermax - iter)/itermax;
    
    % update position
    
    for i = 1:nPop
        mbest = mean(pbest);
        fi    = abs(randn(1,D));
        la    = (c1.*fi.*pbest(i,:)+c2*(1-fi).*gbest./(c1+c2));
        G     = abs(randn(1,D));
        A     = beta.*abs(mbest-pos(i,:)).*log(1./G);
        
        if rand > 0.5
            Xnew = la + A;
        else
            Xnew = la - A;
        end
        % check bound
        Xnew = max(Xnew,lb);
        Xnew = min(Xnew,ub);
        fnew = func1(Xnew);
        
        fx(i) = func1(pos(i,:));
        
        % Update pbest
        
        if fnew < fx(i)
            pbest(i,:) = Xnew;
            pos(i,:)   = Xnew;
            fx(i)      = fnew;
        end
    end
    
% gbest

[fmin,find] = min(fx);

if fmin < Fgbest
    Fgbest = fmin;
    gbest  = pos(find,:);
end

iter = iter+1;
    
%     BestFx = [];
%     BestX  = [];
    
    [optval,optind] = min(fx);
    BestFx(iter)    = optval;
    BestX(iter,:)   = pos(optind,:);
    
    disp(['Iteration' num2str(iter)...
        ': Best Cost = ' num2str(BestFx(iter))]);
    
    % Plot the best
    
    plot(BestFx,'LineWidth',2);
    xlabel('Iteration number');
    ylabel('Fitness value');
    title('Convergence vs Iteration');
    grid on
    
%      stopIter = 10;
%     if iter > stopIter
%         deltaGbestValue = abs(gbest(2:end) - gbest(1:end-1));
%         deltaGbestValueLog = log(deltaGbestValue(deltaGbestValue~=0));
%         if length(deltaGbestValueLog) > 10
%             stopValue = sum(log(gbest(end-9:end)));
%             sumDelta = sum(deltaGbestValueLog(end-9:end));
%         else
%             stopValue = sum(log(gbest));
%             sumDelta = sum(deltaGbestValueLog);
%         end
%         if sumDelta < stopValue/2
%             break;
%         end
%     end
end

study_time = toc;

    

