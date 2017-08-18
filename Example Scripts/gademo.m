% This a script illustrating the use of the Genetic Algorithm (GA)
% optimizer provided with "Engineering Design via Surrogate Modelling".
% As an example it uses the 10-variable Rastrigin function, a classic test
% function, whose main feature is that it has a very large number of local
% optima.
%

%
% Copyright 2009 A Sobester
%
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.

addpath('../Exploring and Exploiting a Surrogate');
addpath('../Example Problems');

NVARS = 10;

global EVALUATIONS
EVALUATIONS = 0;

lb = -5.12*ones(NVARS,1);
ub =  5.12*ones(NVARS,1);

A = [];b = [];Aeq = [];beq = [];

GAOptions = gaoptimset('PopulationSize',2000,...
                       'Generations',35);

[Xstart,bestobj,history,eval_count]=...
    ga(@rastrigin,NVARS,A,b,Aeq,beq,lb,ub,[],GAOptions);

disp('The best optimum found:')
disp(Xstart)

disp('...with an objective value of:')
disp(bestobj)

plot(history)
title('Objective function history')
xlabel('Generation')
ylabel('Objective')