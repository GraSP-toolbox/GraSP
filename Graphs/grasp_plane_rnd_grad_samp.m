%Constructs a complete graph where nodes are non uniformly sampled on the
%2D plane and weighted by a Gaussian kernel.
%
%   graph = GRASP_PLANE_RND_GRAD_SAMP(N) returns a graph constructed from N
%   nodes non uniformly sampled on the plane, using a gradient. The
%   adjacency matrix is weighted by a Gaussian kernel of the distance with
%   parameter sigma = 1. Highest density is in the center.
%
%   GRASP_PLANE_RND(..., options) optional parameters:
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

function graph = grasp_plane_rnd_grad_samp(N, varargin)
    %% Parameters
    default_param = struct(...
        'directed', false,...
        'sigma', 1);
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);

    %% Intializations
    graph = grasp_struct;
    
    %% Random nodes on the plane
    c = 0.3;
    f = @(x) c * x .^ 3 + (1 - c) * x;
    z = (f(2 * rand(N, 1) - 1) * 5) .* exp(1i .* rand(N, 1) * 2 * pi);
    graph.layout = [real(z) imag(z)] + 5;
    
    %% Distances
    graph.distances = grasp_distances_layout(graph);
    
    %% Compute the Gaussian kernel.
    graph.A = grasp_adjacency_gaussian(graph, options.sigma);
    
    %% Symetrise ?
    if ~options.directed
        graph.A = max(graph.A, graph.A');
    end
end