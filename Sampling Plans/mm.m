function Mmplan = mm(X1,X2,p)
% Given two samplig plans chooses the one with the better space-filling
% properties (as per the Morris-Mitchell criterion).
%
% Inputs:
%       X1, X2 - the two sampling plans
%       p - the distance metric to be used (p=1 rectangular - default, p=2
%       Euclidean)
%
% Outputs:
%       Mmplan - if Mmplan = 0, identical plans or equally space-filling,
%                if Mmplan = 1, X1 is more space-filling, if Mmplan = 2, X2
%                is more space-filling.
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

if ~exist('p','var')
    p = 1;
end

if sortrows(X1)==sortrows(X2)
    % If the two matrices contain the same points
    Mmplan = 0;
else
    % Calculate the distance and multiplicity arrays
    [J1,d1] = jd(X1,p); m1 = length(d1);
    [J2,d2] = jd(X2,p); m2 = length(d2);
    
    % Blend the distance and multiplicity arrays together for comparison 
    % according to the Definition 1.4.3B. Note the different signs - we are
    % maximizing the d's and minimizing the J's.
    V1(1:2:2*m1-1) = d1;
    V1(2:2:2*m1)   = -J1;
    
    V2(1:2:2*m2-1) = d2;
    V2(2:2:2*m2)   = -J2;
    
    % The longer vector can be trimmed down to the length of the shorter one
    m = min(m1,m2);
    V1 = V1(1:m); V2 = V2(1:m);

    % Generate vector c such that c(i)=1 if V1(i)>V2(i), c(i)=2 if 
    % V1(i)<V2(i) and c(i)=0 otherwise
    c = (V1 > V2) + 2*(V1 < V2);
    
    % If the plans are not identical but have the same space-filling 
    % properties
    if sum(c)==0
        Mmplan = 0;
    else
        % The more space-filling design (mmplan)
        % is the first non-zero element of c
        i = 1;
        while c(i)==0
            i = i + 1;
        end
        Mmplan = c(i);
    end
end