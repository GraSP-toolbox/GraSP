%Constructs the standard Laplacian of a non directed graph. If the graph is
%directed, with non-negative weights, use Fan Chung combinatorial Laplacian
%[Laplacians and the Cheeger Inequality for Directed Graphs].
%
%   L = GRASP_LAPLACIAN_STANDARD(graph) returns the Laplacian of the
%       graph provided using its adjacency matrix.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018-2019)
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

function L = grasp_laplacian_standard(graph)
    %% Check for self-loops
    if sum(abs(diag(graph.A))) > 0
        error('GraSP:LaplacianStandard:UnsupportedSelfLoops', 'Error: unsupported graph with self-loops!');
    end
    
    %% Compute L
    if grasp_is_directed(graph)
        %% Check that we can convert the adjacency matrix into a random walk matrix P
        if sum(graph.A(:) < 0) > 0
            error('GraSP:LaplacianStandard:UnsupportedNegativeWeights', 'Error: directed graph with negative weights!');
        end
        
        %% Directed graph Laplacian deriving from F. Chung definition
        Dout = grasp_degrees(graph);
        P = Dout ^ -1 * graph.A;
        [pi, ~] = eigs(P', 1, 1); % Closest left eigenvector of P to eigenvalue 1 (stationary distribution of the RW)
        Pi = diag(pi);
        L = Pi - (Pi * P + P' * Pi) / 2;
    else
        %% Classic Combinatorial Laplacian for undirected graph
        L = grasp_degrees(graph) - grasp_adjacency_matrix(graph);
    end
end