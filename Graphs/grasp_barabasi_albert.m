%Constructs a graph from the Barabási–Albert model. The graph is non
%directed.
%
%   graph = GRASP_BARABASI_ALBERT(N, m0) constructs a Barabási–Albert graph
%   with N nodes, and initial degree m0.
%
%   graph = GRASP_BARABASI_ALBERT(N, starting_graph) constructs a
%   Barabási–Albert graph with N nodes, and initialized with
%   starting_graph, and with initial degree equal to the number of nodes of
%   starting_graph.
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

function graph = grasp_barabasi_albert(N, arg2)
    %% Initialization
    graph = grasp_struct;
    graph.A(N, N) = 0;
    
    %% Parameters
    if isstruct(arg2)
        starting_graph = arg2;
        m0 = grasp_nb_nodes(starting_graph);
    else
        m0 = arg2;
        starting_graph = grasp_erdos_renyi(m0, 0.5);
    end
    
    %% Initial graph
    graph.A(1:m0, 1:m0) = starting_graph.A;
    graph.A = (graph.A + graph.A') > 0;
    graph.A = (graph.A + 1) - 1;
    
    %% Intial probability distribution
    cumul = (1:m0) / m0;
    
    %% Add nodes
    for i = (m0 + 1):N
        %% Draw neighbors
        new_neighbors = zeros(1, m0);
        for j = 1:m0
            new_neighbors(j) = 0;
            while new_neighbors(j) == 0 || sum(new_neighbors(1:(j-1)) == new_neighbors(j)) > 0
                r = rand(1);
                new_neighbors(j) = 1;
                while cumul(new_neighbors(j)) <= r
                    new_neighbors(j) = new_neighbors(j) + 1;
                end
            end
        end
        
        %% Add neighbors
        for j = new_neighbors
            graph.A(j, i) = 1;
            graph.A(i, j) = 1;
        end
        
        %% Update probability distribution
        node_cumul_degrees = (blkdiag(tril(ones(i)), zeros(N - i)) * sum(graph.A, 2))';
        cumul = node_cumul_degrees / sum(sum(graph.A));
    end
    
    %% Sparsify
    graph.A = sparse(graph.A);
    
    %% Layout
    graph.layout = grasp_layout(graph);
end