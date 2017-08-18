% Function to be optimized
ObjFun = 'rastrigin';

% Number of design variables
nvars = 10;

% Variable ranges
LB = -5*ones(1,nvars);
UB =  5*ones(1,nvars);

GAOptions.PopulationSize = 3000;
GAOptions.Generations = 20; 
GAOptions.TSize = 6;
GAOptions.Direction = -1; 
GAOptions.Encoding = 12;
GAOptions.Save = 1;
%GAOptions.POP
%GAOptions.Target

[bestvar,bestobj,history,eval_count,ev_hist] = ...
    ga(ObjFun,nvars,[],[],[],[],LB,UB,[],GAOptions);
