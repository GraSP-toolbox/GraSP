%sConstruct a directed torus.
%
%   graph = GRASP_DIRECTED_TORUS(N) constructs an directed torus of NxN,
%   with node (i,j) having id 1 + N * i + j.
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

function graph = grasp_directed_torus(N, M)
    %% Size of the new graph
    if nargin == 1
        M = N;
    end
    NA = N * M;
    
    %% Initialization
    graph = grasp_struct;
    graph.A(NA, NA) = 0;
    graph.layout(NA, 2) = 0;
    
    %% Grid ((0,0) => top left)
    for i = 0:(N - 1)
        for j = 0:(M - 1)
            cur = 1 + M * i + j;
            below = 1 + M * (i + 1) + j;
            right = 1 + M * i + (j + 1);
            if i < N - 1
                graph.A(cur, below) = 1;
            else
                graph.A(cur, 1 + j) = 1;
            end
            if j < M - 1
                graph.A(cur, right) = 1;
            else
                graph.A(cur, 1 + M * i) = 1;
            end
            
            graph.layout(cur, 1) = 10 * (i + 1/2) / N;
            graph.layout(cur, 2) = 10 * (j + 1/2) / M;
        end
    end 
    
    %% Sparsify
    graph.A = sparse(graph.A);
end