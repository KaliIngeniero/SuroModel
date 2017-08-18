% This script illustrates the use of tileplot.m and screeningplot.m
%
% Baseline values
% S_w    -Wing area (ft^2)                              174
% W_fw   -Weight of fuel in the wing (lb)               252
% A      -aspect ratio                                 7.52
% Lambda -quarter-chord sweep (deg)                       0
% q      -dynamic pressure at cruise (lb/ft^2)           34
% lambda -taper ratio                                 0.672
% tc     -aerofoil thickness to chord ratio            0.12
% N_z    -ultimate load factor (1.5x limit load factor) 3.8
% W_dg   -flight design gross weight (lb)             ~2000
% W_p    -paint weight (lb/ft^2)                     ~0.064
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

addpath('../Example Problems');
addpath('../Sampling Plans');
addpath('../Advanced Concepts');

baseline = [174,252,7.52,0 ,34,0.672,0.12,3.8,2000,0.064];
range    = [150,220,6   ,-10,16,0.5  ,0.08,2.5,1700,0.025;
            200,300,10  , 10,45,1    ,0.18,6  ,2500,0.08];
        
labels   = {'\bf S_w','\bf W_{fw}','\bf A','\bf \Lambda','\bf q',...
    '\bf \lambda','\bf tc','\bf N_z','\bf W_{dg}','\bf W_p'};

k = 10;

% Plot the weight function
tileplot(baseline, range, labels, @liftsurfw, 30, 177, 320, 1)

% Screening
r = 5; xi = 1; p = 8;

X = screeningplan(k, p, xi, r);

screeningplot(X, @liftsurfw, range, xi, p, labels);