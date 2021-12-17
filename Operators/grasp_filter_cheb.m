%Defines and returns a Chebyshev polynomial approximation of a given 
%function.
%
%   filter = GRASP_FILTER_CHEB(fun, order) computes the hebyshev polynomial
%       approximation of fun of the first kind of degree order on the
%       interval [0 2].
%
%   GRASP_FILTER_CHEB(..., options) optional parameters:
%
%   options.interval: interval on which to compute the polynomial (default:
%       [0 2]).
%
%   options.cheb_kind: first (1, default), or second (2) of Chebyshev
%       polynomial.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2019)
% 
% benjamin.girault@usc.edu
% 
% This software is a computer program whose purpose is to provide a Matlab
% / Octave toolbox for handling and displaying graph signals.
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.

function filter = grasp_filter_cheb(fun, order, varargin)
    %% Parameters
    default_param = struct(...
        'interval', [0 2],...
        'cheb_kind', 1,...
        'chebfun_splitting', 0);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Filter
    grasp_start_opt_3rd_party('chebfun');
    
    filter.type = 'chebpoly';
    filter.data.coeffs = chebcoeffs(chebfun(fun, options.interval, 'splitting', options.chebfun_splitting), order + 1, 'kind', options.cheb_kind);
    filter.data.interval = options.interval;
    filter.data.cheb_kind = options.cheb_kind;
end