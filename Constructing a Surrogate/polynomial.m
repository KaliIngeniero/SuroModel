function [BestOrder,Coeff,MU] = polynomial(X,Y)
% Fits a one-variable polynomial to one-dimensional data
%
% Inputs:
%       X,Y - training data vectors
%
% Outputs:
%       BestOrder - the order of the polynomial, estimated using
%                   cross-validation
%       Coeff - the coefficients of the polynomial
%       MU - normalisation vector
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


PlotPoints = 100;

OrderLabel =...
    {'Linear',...
    'Quadratic',...
    'Cubic',...
    'Quartic',...
    'Fifth order',...
    'Sixth order',...
    'Seventh order',...
    'Eighth order',...
    'Ninth order',...
    'Tenth order',...
    'Eleventh order',...
    'Twelveth order',...
    'Thirteenth order',...
    'Fourteenth order',...
    'Fifteenth order'};

% The crossvalidation will split the data into this many subsets
% (this may be changed if required)
q = 5;

% This is the highest order polynomial that will be considered
MaxOrder = 15;

n = length(X);

% X split into q randomly selected subsets
XS = randperm(n);

FullXS = XS;

% The beginnings of the subsets...
From = (1:round(n/q):n-1);
To = zeros(size(From));

% ...and their ends
for i=1:q-1
    To(i) = From(i+1)-1;
end

To(q) = n;

CrossVal = zeros(1,MaxOrder);

% Cycling through the possible orders
for Order = 1:MaxOrder
    
    fprintf('Now considering the polynomial of order %d...\n',Order);
    
    CrossVal(Order) = 0;

    % Model fitting to subsets of the data
    for j = 1:q
        Removed = XS(From(j):To(j));
        XS(From(j):To(j)) = [];
       
        [P,S,MU] = polyfit(X(XS),Y(XS),Order);
        
        CrossVal(Order) = CrossVal(Order) +...
            sum((Y(Removed) - polyval(P,X(Removed),S,MU)).^2)...
            /length(Removed);
        XS = FullXS;
    end
    
end

[MinCV, BestOrder] = min(CrossVal);


[Coeff,S,MU] = polyfit(X,Y,BestOrder);


% **** VISUALISING THE TRAINING DATA AND THE PREDICTOR *****

XPrediction = (min(X):(max(X)-min(X))/(PlotPoints-1):max(X));

[YPrediction,DELTA] = polyval(Coeff,XPrediction,S,MU);

plot(X,Y,'o');
hold on
plot(XPrediction,YPrediction,'LineWidth',2);
% Optional error plot
% plot(x,Y+DELTA,'g');
% plot(x,Y-DELTA,'g');
title(strcat(OrderLabel(BestOrder),...
    ' polynomial fitted to the training data'));
xlabel('x');
ylabel('y');
