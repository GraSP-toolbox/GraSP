%Import a graph structure from CSVs.
%
%   graph = GRASP_IMPORTCSV(nodes_file, edges_file) reads the two files
%   nodes_file and edges_file and outputs the graph. This function is the
%   inverse of GRASP_EXPORTCSV.
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

function graph = grasp_importcsv(nodes_file, edges_file)
    graph = grasp_struct;
    
    %% Nodes
    nodes = csvread(nodes_file, 1);
    graph.layout = 10 * nodes(nodes(:, 3), 1:2);
    N = size(graph.layout, 1);
    graph.A = zeros(N);
    
    %% Edges
    edges = csvread(edges_file, 1);
    for i = 1:size(edges, 1)
        graph.A(edges(i, 1), edges(i, 2)) = edges(i, 3);
    end
end