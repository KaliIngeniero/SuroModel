function tileplot(Baseline, Range, Labels, Objhandle,...
                                   Mesh, Lower, Upper, Cont)
%
% Generates the (k-1)x(k-1) tile plot of a k-variable function
%
% Inputs:
%		Baseline - 1xk vector of baseline values assigned to each variable on
%			      a tile where they are not active
%		Range - 2x4 matrix of minimum and maximum values for each
%                         variable
%		Labels - cell array containing the names of the variables
% 		Objhandle - name of the objective function
%		Mesh - the objective function will be sampled on a mesh x mesh full
%              		factorial plan on each tile
% 		Lower/Upper- minimum/maximum value of the objective function - this is 
%              			required to adjust the colour range of each tile with
%              			respect to the full colour range of the function (if not
%              			known, set to [] and the function will estimate it).
% 		Cont- if Cont = 1 contour lines are plotted and the spaces between
%              	        them are filled with a colour determined by the function
%              	        value. Otherwise a colour-shaded plot is generated.
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



% Greying out the end of the colormap - this will be one of the colours
% used to represent infeasible regions. Only available for colour-shaded
% plots (Cont~=1).

if Cont~=1
    Map = colormap;
    Map(end,1:3) = [.5 .5 .5];
    colormap(Map);
end


% Number of variables
k = length(Baseline);

tile_width = 1/k;
tile_height = 1/k;

% The inter-tile gap as a fraction of the width of the full plot
gap = 0.005;

% Data limits for colour axis control
lo = inf;
hi = -inf;

for i=2:k
    for j=1:i-1
        
        rx = Range(1,j):(Range(2,j)-Range(1,j))/(Mesh-1):Range(2,j);
        ry = Range(1,i):(Range(2,i)-Range(1,i))/(Mesh-1):Range(2,i);
        [X,Y] = meshgrid(rx,ry);
        
        
        % Computing tile i,j
        M = zeros(size(X));
        for l = 1:size(X,1)
            for ll = 1:size(X,2)
                design = Baseline;
                design(j) = X(l,ll);
                design(i) = Y(l,ll);
                M(l,ll) = feval(Objhandle,design);
                if M(l,ll)>hi, hi=M(l,ll); end
                if M(l,ll)<lo, lo=M(l,ll); end                
            end
        end
        
        % Plotting the current tile
        axes;
        hold on
        set(gca,'Position',[tile_width/2+(j-1)*tile_width,...
            tile_height/2+(k-i)*tile_height,...
            tile_width-gap,tile_height-gap])
        
        
        if Cont == 1
            contourf(X,Y,M)
            % If you don't want the contourlines to show, use:
            % contourf(X,Y,M,'LineStyle','none')
            % instead of the line above
        else
            pcolor(X,Y,M), shading interp, axis tight
        end
        
     
        
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        set(gca,'FontSize',12)
        if ~isempty(Lower) && ~isempty(Upper)
            caxis([Lower Upper])
        else
            caxis([lo hi]);
            Upper = hi;
            Lower = lo;
        end
        if j==1
            ylabel(Labels(i));
        end
        if i==k
            xlabel(Labels(j));
        end
    end
end

fprintf('The extrema of the evaluated objective function values are:\n');
fprintf('Minimum: %4.3f\n',lo);
fprintf('Maximum: %4.3f\n',hi);

% If the user-specified objective function value range (which is needed to
% specify the colour range for each tile) is outside the actual objective
% function range by more than 5% either side...
if abs(Lower-lo)>abs(0.05*lo) || abs(Upper-hi)>abs(0.05*hi)
    fprintf('You may wish to re-run the function with these values\n');
    fprintf('as the last two arguments.\n');
end
    