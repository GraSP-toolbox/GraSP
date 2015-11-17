%Export a graph structure in CSV files for LaTeX.
%
%   GRASP_EXPORTCSV(graph, nodes_file, edges_file) creates the two files
%   nodes_file and edges_file that can be used by LaTeX package stemgraph.
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

function grasp_exportcsv(graph, nodes_file, edges_file)
    N = grasp_nb_nodes(graph);
    fileID = fopen(nodes_file,'w');
    fprintf(fileID,'X,Y,Name\n');
    for i = 1:N
        fprintf(fileID, '%d,%d,%d\n', graph.layout(i, 1) / 10, graph.layout(i, 2) / 10, i);
    end
    fclose(fileID);

    fileID = fopen(edges_file, 'w');
    fprintf(fileID,'n1,n2,w12,oriented\n');
    for i = 1:N
        for j = 1:N
            directed = full(graph.A(i,j) ~= graph.A(j,i));
            if ~directed && i > j
                continue;
            end
            if (abs(graph.A(i, j)) > 0)
                fprintf(fileID, '%d,%d,%e,%d\n', i, j, full(graph.A(i, j)), directed);
            end
        end
    end
    fclose(fileID);
end