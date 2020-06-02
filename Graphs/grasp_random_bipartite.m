%Constructs a random bipartite graph.
%Attention: the resulting graph may not be connected!
%The y coordinates of the vertices is random, while the x coordinate is
%either 3 or 7 depending on the which set the vertex belongs to.
%
%   graph = GRASP_RANDOM_BIPARTITE() constructs a random bipartite graph
%   with 50 vertices in each set.
%
%   GRASP_RANDOM_BIPARTITE(N) constructs a random bipartite graph with N
%   vertices in each set.
%
%   GRASP_RANDOM_BIPARTITE(N, M) constructs a random bipartite graph with N
%   vertices in the first set (left on the figure) and M in the second
%   (right on the figure).
%
%   GRASP_RANDOM_BIPARTITE(..., options) optional parameters:
%
%   options.p: probability to draw an edge (default: 2 * log(k) / k with
%   k = max(N, M)).
%
% Authors:
%  - Sunil K Narang <kumarsun@usc.edu>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Sunil K Narang, University of Sourthern California, Los
% Angeles, California, USA (2011)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016)
% 
% kumarsun@usc.edu
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

function graph = grasp_random_bipartite(varargin)
    %% Parameters
    options = struct(...
        'N', 50,...
        'M', 50,...
        'p', -1);
    if nargin == 1
        if isnumeric(varargin{1})
            options.N = varargin{1};
        else
            options = grasp_merge_structs(options, varargin{1});
        end
    elseif nargin == 2
        if isnumeric(varargin{2})
            options.N = varargin{1};
            options.M = varargin{2};
        elseif isnumeric(varargin{1})
            options = grasp_merge_structs(options, varargin{2});
            options.N = varargin{1};
            options.M = varargin{1};
        else
            options = grasp_merge_structs(options, cell2struct(varargin(2:2:end), varargin(1:2:end), 2));
        end
    elseif nargin == 3
        options = grasp_merge_structs(options, varargin{3});
        options.N = varargin{1};
        options.M = varargin{2};
    elseif nargin > 3
        options = grasp_merge_structs(options, cell2struct(varargin(4:2:end), varargin(3:2:end), 2));
        options.N = varargin{1};
        options.M = varargin{2};
    end
    
    if options.p == -1
        max_set = max(options.N, options.M);
        options.p = 2 * log(max_set) / (max_set);
    end
    
    %% Graph construction
    graph = grasp_struct;
    graph.A = zeros(options.N + options.M);
    graph.A(1:options.N, (options.N + 1):(options.N + options.M)) = rand(options.N, options.M) < options.p;
    graph.A = graph.A + graph.A';
    
    %% Layout
    graph.layout(1:(options.N + options.M), 2) = 10 * rand(options.N + options.M, 1);
    graph.layout(1:options.N, 1) = 3;
    graph.layout(options.N + (1:options.M), 1) = 7;
end