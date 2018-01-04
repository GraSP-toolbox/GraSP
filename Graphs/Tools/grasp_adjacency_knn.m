%Constructs an adjacency matrix using the k nearest neighbors, assuming
%that the adjacency matrix is a similarity matrix (bigger weights for closer
%nodes).
%
%   A = GRASP_ADJACENCY_KNN(graph, k) returns an adjacency matrix A
%   constructed from the adjacency matrix using KNN.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Copyright Benjamin Girault, University of Southern California, USA
% (2017).
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

function A = grasp_adjacency_knn(graph, k)
    %% Intializations
    N = grasp_nb_nodes(graph);
    directed = grasp_is_directed(graph);
    
    %% Compute the k nearest neighbors
    I(2 * k * N) = 0;
    J(2 * k * N) = 0;
    W(2 * k * N) = 0;
    nb_entries = 0;
    vertex_marking = zeros(N);
    for i = 1:N
        mask = [1:(i-1) (i+1):N];
        [~, IX] = sort(graph.A(i, mask), 'descend');
        IX = mask(IX);
        for j = IX(1:k)
            if vertex_marking(i, j) == 0
                I(nb_entries + 1) = i;
                J(nb_entries + 1) = j;
                W(nb_entries + 1) = graph.A(i, j);
                vertex_marking(i, j) = 1;
                nb_entries = nb_entries + 1;
            end
            if ~directed && vertex_marking(j, i) == 0
                I(nb_entries + 1) = j;
                J(nb_entries + 1) = i;
                W(nb_entries + 1) = graph.A(j, i);
                vertex_marking(j, i) = 1;
                nb_entries = nb_entries + 1;
            end
        end
    end
    A = sparse(I(1:nb_entries), J(1:nb_entries), W(1:nb_entries), N, N);
end