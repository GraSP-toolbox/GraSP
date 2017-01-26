%Constructs a graph of k the nearest neighbors where nodes are sampled on
%the 2D torus.
%
%   graph = GRASP_TORUS_KNN(N, k) returns a graph constructed from nodes
%   uniformly sampled on the 2D torus and connected using the k nearest
%   neighbors method. The adjacency matrix is weighted by the Gaussian
%   kernel of the distance.
%
%   GRASP_TORUS_KNN(..., options) optional parameters:
%
%   options.directed: true if graph should stay directed (default: false)
%   options.sigma: sigma^2 in the Gaussian kernel (default: 1)
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% 
% benjamin.girault@ens-lyon.fr
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

function graph = grasp_torus_knn(N, k, varargin)
    %% Parameters
    default_param = struct(...
        'directed', false,...
        'sigma', 1);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);

    %% Intializations
    graph = grasp_struct;
    
    %% Random nodes on the plane
    graph.layout = 10 * rand(N, 2);
    
    %% Distances
    [dx1, dx2] = meshgrid(graph.layout(:, 1), graph.layout(:, 1));
    dx = min(dx2 - dx1, min(dx1 + 10 - dx2, 10 - dx1 + dx2));
    [dy1, dy2] = meshgrid(graph.layout(:, 2), graph.layout(:, 2));
    dy = min(dy2 - dy1, min(dy1 + 10 - dy2, 10 - dy1 + dy2));
    graph.distances = sqrt(dx .^ 2 + dy .^ 2);
    
    %% Compute the k nearest neighbors
    graph.A = grasp_adjacency_gaussian(graph, options.sigma);
    graph.A = grasp_adjacency_knn(graph, k);
    
    %% Symetrise ?
    if ~options.directed
        graph.A = max(graph.A, graph.A');
    end
end