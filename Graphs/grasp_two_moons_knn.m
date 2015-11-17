%Constructs a two moons graph with the k nearest neighbors.
%
%   graph = GRASP_TWO_MOONS_KNN(N, sigma_noise, k) returns a graph
%   constructed from N nodes sampled around two moons. Uses a gaussian
%   noise of standard deviation sigma_noise. The edges are drawn with the k
%   nearest neighbors approach.
%
%   GRASP_TWO_MOONS_KNN(..., options) optional parameters:
%
%   options.sigma: sigma^2 in the Gaussian kernel
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

function graph = grasp_two_moons_knn(N, sigma_noise, k, varargin)
    %% Parameters
    default_param = struct(...
        'sigma', 1);
    if nargin == 3
        options = struct;
    elseif nargin > 4
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);

    %% Intializations
    graph = grasp_struct;
    graph.A(N, N) = 0;
    
    %% Two moons
    N1 = round(2 * N / 3);
    t1 = 2 * rand(N1, 1);
    t2 = 2 * rand(N - N1, 1) + 2;
    t = [t1 ; t2];
    graph.class(1:N1) = 1;
    graph.class((N1 + 1):N) = 2;
    graph.layout = moons(t);
    
    %% Noise
    graph.layout = graph.layout + normrnd(0, sigma_noise, size(graph.layout));
    
    %% Distances
    graph.distances = grasp_distances_layout(graph);
    
    %% Compute the k nearest neighbors
    graph.A = grasp_adjacency_gaussian(graph, options.sigma);
    graph.A = grasp_adjacency_knn(graph, k);
    graph.A = max(graph.A, graph.A');
end

function out = moons(t)
    p = 3;
    
    %% Upper moon
    x = (4 + p * (t - 3)) .* (t > 2);
    y = (6.1 - p / 2 * (t - 3) .^ 2) .* (t > 2);
    
    %% Lower moon
    x = x + (6 + p * (t - 1)) .* (t <= 2);
    y = y + (3.9 + p / 2 * (t - 1) .^ 2) .* (t <= 2);
    
    out = [x y];
end