function rbf
% Estimates the parameters of a Radial Basis Function model (all other
% information about the model is expected to exist in the ModelInfo global
% variable).
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

global ModelInfo

if ModelInfo.Code < 4
	% Fixed basis function, only w needs estimating
	rbfestimw
else
	% Basis function also requires a sigma,
	% this needs to be estimated first
	MI = ModelInfo;

	% Direct search between 10^-2 and 10^2
    Sigmas = logspace(-2,2,30);

	% Number of cross-validation subsets
	if size(MI.X,1) < 6
		q = 2;
        elseif size(MI.X,1) < 15
            q = 3;
        elseif size(MI.X,1) < 50
            q = 5;
        else
            q = 10;
	end
		
	% Number of sample points
	n = size(ModelInfo.X,1);

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

	CrossVal = zeros(1,length(Sigmas));

	% Cycling through the possible values of Sigma
	for Sigindex = 1:length(Sigmas)
    
        fprintf('Computing cross-validation metric for Sigma=%6.4f...\n',...
		Sigmas(Sigindex));
    
		CrossVal(Sigindex) = 0;

		% Model fitting to subsets of the data
    	for j = 1:q
        	Removed = XS(From(j):To(j));
    		XS(From(j):To(j)) = [];

  
            ModelInfo.X = MI.X(XS,:);
        	ModelInfo.y = MI.y(XS);
    		ModelInfo.Sigma = Sigmas(Sigindex);
        
        	% Sigma and subset chosen, now estimate w
    		rbfestimw
        
            if isempty(ModelInfo.Weights)
                CrossVal(Sigindex) = 1e20;
                XS = FullXS;
                break
            end
            
            % Compute vector of predictions at the removed sites
    		Pr = zeros(length(Removed),1);
    		for jj = 1:length(Removed)
    			Pr(jj) = predrbf(MI.X(Removed(jj),:));
            end
        
        
    		CrossVal(Sigindex) = CrossVal(Sigindex) +...
        		sum((MI.y(Removed) - Pr).^2)/length(Removed);
        
        
    		XS = FullXS;
        end


        % Now attempt Cholesky on the full set, in case the subsets could
        % be fitted correctly, but the complete X could not
        ModelInfo = MI;
        ModelInfo.Sigma = Sigmas(Sigindex);
        rbfestimw
        if isempty(ModelInfo.Weights)
            CrossVal(Sigindex) = 1e20;
            display('Failed to fit complete sample data.')
        end



    end

	[MinCV, BestSig] = min(CrossVal);



	ModelInfo = MI;
	fprintf('Selected sigma=%5.4f\n',Sigmas(BestSig))
	ModelInfo.Sigma = Sigmas(BestSig);
	
	rbfestimw

end
	


function rbfestimw
% Estimates the basis function weights w if there is no other parameter
% (sigma), or this is known already and included in ModelInfo

global ModelInfo


% This should have the following fields defined:
% Training data: ModelInfo.X
%                ModelInfo.y
% Basis function type:
%                ModelInfo.Code
% If basis function type > 3
% additional basis function parameter
%                ModelInfo.Sigma


if (ModelInfo.Code > 3) && (~isfield(ModelInfo,'Sigma'))
    error('The basis function in ModelInfo.Code requires a Sigma.')
end

%Number of points
n = length(ModelInfo.y);

% Build Gram matrix
d = zeros(n,n);
for i=1:n
   for j=1:i
      d(i,j) = norm(ModelInfo.X(i,:)-ModelInfo.X(j,:),2);      
      d(j,i) = d(i,j);
   end
end

% Constructing the PHI matrix
ModelInfo.Phi = zeros(n,n);
for i = 1:n
   for j = 1:i
      if isfield(ModelInfo,'Sigma')
          ModelInfo.Phi(i,j) = basis(ModelInfo.Code, d(i,j), ModelInfo.Sigma);
      else
          ModelInfo.Phi(i,j) = basis(ModelInfo.Code, d(i,j));
      end
      ModelInfo.Phi(j,i) = ModelInfo.Phi(i,j);
   end
end

% Calculating the weights of the radial basis function surrogate
% Cholesky factorization used if Gaussian or inverse mq required
% LU decomposition otherwise

if ModelInfo.Code==4 ||  ModelInfo.Code==6
    [U,p] = chol(ModelInfo.Phi);
    if p==0
        ModelInfo.Weights = U\(U'\ModelInfo.y);
        ModelInfo.Success = 1;
    else
        display('Cholesky factorization failed.');
        display('Two points may be too close together.');
        ModelInfo.Weights = [];
        ModelInfo.Success = 0; 
    end
else
    ModelInfo.Weights = ModelInfo.Phi\ModelInfo.y;
    if isempty(ModelInfo.Weights)
        ModelInfo.Success = 0;
    else
        ModelInfo.Success = 1;
    end
end
