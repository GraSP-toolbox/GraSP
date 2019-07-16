%Returns a cell array of subgraphs that are maximally connected.
%
%   graphs = GRASP_MAXIMALLY_CONNECTED_SUBGRAPHS(graph) compute the
%       maximally connected subgraphs of graph and return them in the
%       struct array graphs.
%
%   [..., nodes] = GRASP_MAXIMALLY_CONNECTED_SUBGRAPHS(...) also returns
%       the set of selected nodes in the original graph for each subgraph.
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

function [graphs, nodes] = grasp_maximally_connected_subgraphs(graph)
    % Load the matlab-bgl toolbox for 'components'
    grasp_start_opt_3rd_party('MatlabBGL');
    
    % Find the connected component
    [conn_comps, sizes] = components(sparse(graph.A));
    
    % Sort by decreasing size
    [~, IX] = sort(sizes, 'descend');
    
    % Build the graph by decreasing size
    graphs(numel(IX)) = graph;
    nodes = cell(numel(IX), 1);
    for i = 1:numel(IX)
        nodes{i} = find(conn_comps == IX(i));
        grasp_subgraph(graph, nodes{i})
        graphs(i)
        graphs(i) = grasp_subgraph(graph, nodes{i});
    end
    
%     [~, largest_cc] = max(sizes);
%     nodes = find(conn_comps == largest_cc);
%     
%     % Select only the nodes of that component
%     graph = grasp_subgraph(graph, nodes);
end