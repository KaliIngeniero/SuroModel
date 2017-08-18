function screening_plot(X, Objhandle, Range, xi, p, Labels)
% Generates a variable elementary effect screening plot
%
% Inputs:
%       X - screening plan built within a [0,1]^k box (e.g. with
%           screening_plan.m)
%       Objhandle - name of the objective function
%       Range - 2xk matrix (k - number of design variables) of lower bounds
%               (first row) and upper bounds (second row) on each variable.
%       xi- elementery effect step length factor
%       p - number of discreet levels along each dimension
%       Labels - 1xk cell array containing the names of the variables
%
% Copyright 2007 A Sobester
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

k = size(X,2);
r = size(X,1)/(k+1);


for i=1:size(X,1)
    X(i,:) = Range(1,:) + X(i,:).*(Range(2,:)-Range(1,:));
    t(i) = feval(Objhandle,X(i,:));
end

for i=1:r
    for j = (i-1)*(k+1)+1:(i-1)*(k+1)+k
       F(find(X(j,:)-X(j+1,:)~=0),i) = (t(j+1)-t(j))/(xi/(p-1));
    end
end

% Compute statistical measures
for i=1:k
    ssd(i) = std(F(i,:));
    sm(i)  = mean(F(i,:));
end
 
figure, hold on
 
for i=1:k
    text(sm(i),ssd(i),Labels(i),'FontSize',25)
end
 
axis([min(sm) max(sm) min(ssd) max(ssd)]);
xlabel('Sample means')
ylabel('Sample standard deviations')
set(gca,'FontSize',14)   