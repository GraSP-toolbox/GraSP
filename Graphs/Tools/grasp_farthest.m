%Determines the two farthest nodes of a graph.
%
%   [x,y] = GRASP_FARTHEST(graph) returns the two farthest nodes of graph.
%
%   GRASP_FARTHEST(..., 'adjacency_matrix', am) type of the adjacency
%   matrix of graph. Default is 'distance' (graph.A is a distance matrix).
%   Other possible values are:
%     * 'similarity' (the distance used is e^-A)
%     * 'gaussian_kernel' (the kernel is inverted to get the distance)
%     * 'dst_mat' (graph.distance will be used).
%
%   [x,y,distances] = GRASP_FARTHEST(...) also returns the matrix of
%   distances between every couple of points.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018)
% 
% benjamin.girault@ens-lyon.fr
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

function [x, y, distances] = grasp_farthest(graph, varargin)
    %% Parameters
    default_param = struct(...
        'adjacency_matrix', 'distance');
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Load the matlab-bgl toolbox for 'all_shortest_paths'
    grasp_start_opt_3rd_party('MatlabBGL');
    
    %% Edge-distance matrix
    switch options.adjacency_matrix
        case 'distance'
            distances = graph.A;
        case 'similarity'
            distances = exp(-graph.A) - eye(grasp_nb_nodes(graph));
        case 'gaussian_kernel'
            if issparse(graph.A)
                distances = spfun(@(x) sqrt(-log(x)), graph.A);
            else
                distances = arrayfun(@(x) sqrt(-log(x)), graph.A);
            end
        case 'dst_mat'
            distances = graph.distances;
        otherwise
            error('Unrecognized ''adjacency_matrix'' argument!');
    end
    
    %% Matrix of distances in the graph
    distances = all_shortest_paths(sparse(distances));
    
    %% Get the two farthest nodes
    [column_max, column_max_ids] = max(distances);
    [~, y] = max(column_max);
    x = column_max_ids(y);
end