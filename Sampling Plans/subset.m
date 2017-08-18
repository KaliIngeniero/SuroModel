function [Xs,Xr] = subset(X,ns)
% Given a sampling plan, returns a subset of a given size with optimized
% space-filling properties (as per the Morris-Mitchell criterion).
%
% Inputs:
%       X - full sampling plan
%       ns - size of the desired subset
%
% Outputs:
%       Xs - subset with optimized space-filling properties
%       Xr - remainder X\Xs
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

n=size(X,1);

% Norm and quality metric exponent - different values can be used if
% required
p = 1; q = 5;

r = randperm(n);

Xs=X(r(1:ns),:);
Xr=X(r(ns+1:end),:);

for j=1:ns
    orig_crit = mmphi(Xs,q,p);
    orig_point = Xs(j,:);
          
    % We look for the best point to substitute the current one with
    bestsub = 1;
    bestsubcrit = Inf;
    
    for i=1:n-ns
        % We replace the current, jth point with each of the remaining
        % points, one by one
        Xs(j,:)=Xr(i,:);
        crit = mmphi(Xs,q,p);
        
        if crit < bestsubcrit
            bestsubcrit = crit;
            bestsub = i;
        end

    end

    if bestsubcrit < orig_crit
        Xs(j,:) = Xr(bestsub,:);
        Xr(bestsub,:)=orig_point;
    else
        Xs(j,:) = orig_point;
    end
            
end
