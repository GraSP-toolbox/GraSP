%Computes a layout of a graph using its spectral properties.
%Uses the random walk Laplacian eigenvectors to build the embedding.
%
%   layout = GRASP_LAYOUT_SPECTRAL(graph) given a graph, computes a layout
%       (i.e. a n-by-2 matrix) of the nodes.
%
%   GRASP_LAYOUT_SPECTRAL(..., options) optional parameters:
%
%   options.normalize_coordinates: whether or not the coordinates returned 
%       need to be normalized to the square [0,10] (default: true).
%
%   options.dimension: dimension of the layout (default: 2 for 2D).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2017-2019)
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

function layout = grasp_layout_spectral(graph, varargin)
    %% Parameters
    default_param = struct(...
        'normalize_coordinates', true,...
        'dimension', 2);
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Main code
    N = grasp_nb_nodes(graph);
    [layout, ~] = eigs(graph.A, grasp_degrees(graph), options.dimension + 1);
    layout = layout(:, 2:(options.dimension + 1));
    if nargin == 1 || options.normalize_coordinates
        layout = 10 * (layout - ones(N, 1) * min(layout)) ./ (ones(N, 1) * (max(layout) - min(layout)));
    end
end