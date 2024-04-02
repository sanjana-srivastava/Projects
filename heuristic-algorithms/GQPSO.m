function [x,fval,exitFlag,output]=GQPSO(fun,nvars,lb,ub)
%% QPSO Parameters
D=16;                      %Number of dimensions
n=100;                     %Population size                  
itermax=100;               %Maximum number of iterations
c1=1;                      %constants
c2=1;
w1=0.5;
w2=1;

%% Generating initial population
for i=1:n
    for j=1:D
        pos(i,j)=lb(j)+rand.*(ub(j)-lb(j));   
    end
end

%% Evaluate objective function
fx=fun(pos);

%% Initialize pbest & gbest
pbest=pos;
[fminvalue,ind]=min(fx);
gbest=pos(ind,:);

%%QPSO Main Loop
iter=1;
while iter<itermax
    Fgbest=fminvalue;
    bet=w1+(w2-w1)*(itermax-iter)/itermax;
    for i=1:n
        mbest=mean(pbest);
        fi=abs(randn(1,D));
        la=(c1*fi.*pbest(i,:)+c2*(1-fi).*gbest)./(c1+c2);
        G=abs(randn(1,D));
        A=bet.*abs(mbest-pos(i,:)).*log(1./G);
        if rand>0.5
            Xnew=la+A;
        else
            Xnew=la-A;
        end
        Xnew=max(Xnew,lb);
        Xnew=min(Xnew,ub);
        fnew=fun(Xnew);
        if fnew<fx(i)
            pbest(1,:)=Xnew;
            pos(i,:)=Xnew;
            fx(i)=fnew;
        end
    end
    [fmin,find]=min(fx);
    if fmin<Fgbest
        Fgbest=fmin;
        gbest=pos(find,:);
    end
    iter=iter+1;
    [optval, optind]=min(fx);
    BestFx(iter)=optval;
    BestX(iter,:)=pos(optind,:);
    disp(['Iteration ' num2str(iter) ...
        ': Objective Value = ' num2str(BestFx(iter))]);
    plot(BestFx,'Linewidth',2);
    xlabel('Iteration Number');
    ylabel('FitnessValue')
    title('Convergence vs Iteration')
    grid on
end

    
    
        
         
        