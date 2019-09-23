%Constructs a complete graph where nodes are sampled on the 2D plane and
%weighted by a Gaussian kernel. Nodes will be scattered pseudo-regularly
%on the 2D plane.
%
%   graph = GRASP_PLANE_RND(N, M) returns a graph constructed from M nodes
%       uniformly sampled on the plane, and with one of the closest 
%       vertices iteratively removed until N remain. The adjacency matrix
%       is weighted by a Gaussian kernel of the distance with parameter 
%       sigma = 1.
%
%   GRASP_PLANE_RND(..., options) optional parameters:
%
%   options.directed: true if graph should stay directed (default: false)
%
%   options.sigma: sigma^2 in the Gaussian kernel (default: 1)
%
%   options.display_threshold: edge weight threshold for display purposes
%       (default:0.3)
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

function [graph, graph_orig] = grasp_plane_rnd_pseudo_regular(N, M, varargin)
    %% Parameters
    default_param = struct(...
        'directed', false,...
        'sigma', 1,...
        'display_threshold', 0.3);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);

    %% Intializations
    graph = grasp_plane_rnd(M);
    graph_orig = graph;
    
    %% Trim vertices
    maxi = max(graph.distances(:));
    cur_N = M;
    while cur_N > N
        tmp = graph.distances + 2 * maxi * eye(grasp_nb_nodes(graph));
        [~, i] = min(tmp(:));
        [i, j] = ind2sub(size(graph.distances), i);
        if sum(graph.A(:, i)) < sum(graph.A(:, j))
            i = j;
        end
        graph = grasp_subgraph(graph, [1:(i-1) (i+1):cur_N]);
        cur_N = grasp_nb_nodes(graph);
    end
    
    %% Compute the Gaussian kernel.
    graph.A = grasp_adjacency_gaussian(graph, options.sigma);
    
    %% Symetrise ?
    if ~options.directed
        graph.A = max(graph.A, graph.A');
    end
    
    %% Threshold the weights to keep the number of displayed edges small
    graph.A_layout = grasp_adjacency_thresh(graph, options.display_threshold);
end