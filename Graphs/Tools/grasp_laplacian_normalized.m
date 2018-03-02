%Constructs the normalized Laplacian of the (possibly directed) graph.
%
%Note: The graph provided can be directed, but then all weights must be
%non-negative.
%
%   Gamma = GRASP_LAPLACIAN_NORMALIZED(graph) returns the Laplacian of the
%   graph provided using its adjacency matrix (see [Li & Zhang 2012:
%   Digraph Laplacian and the Degree of Asymmetry])
%
%   [Gamma, Phi] = GRASP_LAPLACIAN_NORMALIZED(...) returns also the
%   stationnary distribution Phi (diagonal matrix)
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Angeles, California, USA (2017-2018)
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

function [Gamma, Phi] = grasp_laplacian_normalized(graph)
    %% Initializations
    N = grasp_nb_nodes(graph);
    D = grasp_degrees(graph);
    P = D^-1 * graph.A;

    %% Checks and construction of Phi    
    if ~grasp_is_directed(graph)
        %% Non directed graph: Phi = D
        Phi = D;
        Dsqrtinv = diag(diag(D) .^ (-1/2));
        Gamma = Dsqrtinv * (D - graph.A) * Dsqrtinv;
        return;
    else
        %% Directed graph: Phi = Froebinius vector
        if sum(sum(graph.A < 0)) > 0
            error('Error: negative weights with a directed graph!');
        end
            
        [V, ~] = eig(P');
        phi = V(:, sum(V >= 0, 1) == N);
        if size(phi, 2) == 0
            phi = V(:, sum(V <= 0, 1) == N);
            phi = -phi;
        end
        Phi = diag(phi);
    end
    
    %% Diplacian
    Phi = full(Phi);
    Gamma = Phi ^ (1/2) * (speye(N) - P) * Phi ^ (-1/2);
    Gamma = sparse(Gamma);
end