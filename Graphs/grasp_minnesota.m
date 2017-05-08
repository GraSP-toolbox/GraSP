%Construct the Minnesota transport graph.
%
%   graph = GRASP_MINNESOTA() returns the Minnesota transport graph of
%   3rdParty/MathBFL.
%
%   GRASP_MINNESOTA(options) optional parameters:
%
%   options.type: specify how to construct the edges, 'gauss_dist' for a
%       Gaussian kernel of the distance (default), 'shortest_path' for an
%       adjacency matrix based on a Gaussian kernel of the shortest path
%       (road) distance, 'road_dist' for the each edge corresponding to a
%       road and weighted by its length, 'road' for 0/1 edges corresponding
%       to roads.
%   options.sigma: sigma^2 in the Gaussian kernel
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016)
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

function graph = grasp_minnesota(varargin)
    %% Parameters
    default_param = struct(...
        'type', 'gauss_dist',...
        'sigma', 1);
    if nargin == 0
        options = struct;
    elseif nargin > 1
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Initializations
    graph = grasp_struct;
    data = load('minnesota');
    % Remove nodes #348 and #349 (disconnected nodes)
    mask = [1:347 350:numel(data.labels)];
    
    %% Basic data
    graph.A = data.A(mask, mask);
    graph.node_names = data.labels(mask);
    graph.layout = data.xy(mask, :);
    
    %% Layout rescaled to the [0;10] square
    m = ones(size(graph.layout, 1), 1) * min(graph.layout);
    M = ones(size(graph.layout, 1), 1) * max(graph.layout);
    graph.layout = 10 * (graph.layout - m) ./ (M - m);
    
    %% Distance matrix from the layout (between any two nodes)
    graph.distances = grasp_distances_layout(graph);
    
    %% Adjacency matrix
    graph.A_layout = graph.A;
    switch options.type
        case 'gauss_dist'
            % nothing yet to do
        case 'shortest_path'
            % changing the distances
            graph.distances = all_shortest_paths(graph.distances .* graph.A);
            graph.distances = (graph.distances + graph.distances') / 2;
        case 'road_dist'
            % changing weights to distances, and we are done
            graph.A = graph.distances .* graph.A;
            return;
        case 'road'
            % nothing more to do
            return;
        otherwise
            error('Unrecognized type of graph! Please read the documentation for valid types.');
    end
    graph.A = grasp_adjacency_gaussian(graph, options.sigma);
end