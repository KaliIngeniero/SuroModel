% Another example of using nested4.m, this time using the
% lifting surface weight example
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
addpath('../Advanced Concepts');

range    = [150,0.08,2.5,1700;
            200,0.18,6  ,2500];
        
labels   = {'\bf S_w','\bf W_{fw}','\bf A','\bf \Lambda','\bf q',...
    '\bf \lambda','\bf tc','\bf N_z','\bf W_{dg}','\bf W_p'};

objhandle = 'liftsurfw4';
        

nested4([3 2 1 4], [10 10], range, labels, objhandle,...
                                   100, 151.17117, 414.20208, 0)