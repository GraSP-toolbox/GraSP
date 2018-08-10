%Constructs an Erdős-Renyi graph.
%
%   graph = GRASP_ERDOS_RENYI(N, p) returns the adjacency matrix A of a
%       Erdős-Renyi graph of N nodes with probability p on the edges. The
%       graph is directed.
%
%   GRASP_ERDOS_RENYI(..., options) optional parameters:
%
%   options.directed: true if graph should stay directed (default: true)
%
%   options.build_layout: true if the graph should include a layout
%       (default: true)
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

function graph = grasp_erdos_renyi(N, p, varargin)
    %% Parameters
    default_param = struct(...
        'directed', true,...
        'build_layout', true);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Graph Realization
    graph = grasp_struct;
    
%     graph.A = rand(N, N) <= p;
    graph.A = rand(N, N);
    if ~options.directed
        % This method is faster than triu then A + A'
        if p <= 0.5
            graph.A = double((graph.A + graph.A') < sqrt(2 * p));
        else
            graph.A = double((graph.A + graph.A') < 2 - sqrt(2 - 2 * p));
        end
    else
        graph.A = double(graph.A <= p);
    end
    graph.A = graph.A - diag(diag(graph.A));
    
    %% Layout?
    if options.build_layout
        graph.layout = grasp_layout(graph);
    end
end