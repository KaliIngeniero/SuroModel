% One variable SVR example (produces figure 2.12)
%
% Copyright 2007 A I J Forrester
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

rand('state',0)
randn('state',0)

% Sample data
X=0:0.05:1;
n=length(X);
X=X';
y=(6.*X-2).^2.*sin(12.*X-4)+randn(n,1);

% Build correlation matrix (set sigma=30 arbitrarily without tuning)
sigma=0.2;
Psi=zeros(n,n);
for i=1:n
    for j=1:n
        Psi(i,j)=exp(-(1/sigma^2)*(X(i)-X(j))^2);
    end
end

% User defined constants
e=1;
C=1e3;
xi=1e-6;

% Construct dual variable optimization problem
% Matrix of correlations
Psi=[Psi -Psi;-Psi Psi];
% Constraint terms
c=[(e*ones(n,1) - y); (e*ones(n,1) + y)]; 
% Lower bound |alpha|>=0
lb=zeros(2*n,1);
% Upper bound |alpha|<=0
ub=(C/n)*ones(2*n,1);
% Start at alpha=[0;0;...;0]
x0=zeros(2*n,1);
% Set sum(alpha^+ - alpha^-)=0
Aeq=[-ones(1,n) ones(1,n)];
beq=0;

% Run quadprog
alpha=quadprog(Psi,c,[],[],Aeq,beq,lb,ub,x0);

% Combine alphas into nx1 vector of SVs
alpha_pm=alpha(1:n)-alpha(n+1:2*n)
% Find indices of SVs
sv_i=find(abs(alpha_pm)>xi);
num_svs=length(sv_i)

% Find SV mid way between 0 and C for mu calculation
[sv_mid,sv_mid_i]=min(abs(abs(alpha_pm)-(C/(2*n))));
% Calculate mu
mu=y(sv_mid_i)-e*sign(alpha_pm(sv_mid_i))-alpha_pm(sv_i)'*Psi(sv_i,sv_mid_i);

% Calculate predictions for plotting
x=[0:0.01:1]; 
for i = 1:length(x)
    % Basis function vaules at point to be predicted
    for j = 1:n
        psi(j,1)=exp(-(1/sigma^2)*(x(i)-X(j))^2);
    end
    % SVR prediction
	pred(i)=mu+alpha_pm'*psi;
end

% Plot prediction
plot(x,pred,'k','LineWidth',2)
hold on
% Plot e-tube
plot(x,pred+e,'k--')
plot(x,pred-e,'k--') 
for i = 1:size(y)
    % Plot sample data
    plot(X(i),y(i),'k+') 
    % Circle support vectors
    if (abs(alpha_pm(i)) > xi)
        plot(X(i),y(i),'ko','markersize',10) 
	end
end
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)
l=legend('prediction','+\epsilon','-\epsilon','sample data','support vectors',2)
set(l,'FontSize',16)
set(gca,'FontSize',16)
axis([0 1 -10 17])
