%Constructs the adjacency matrix of a graph as a Gaussian kernel of its
%distance matrix. If the graph already has an adjacency matrix, only the
%previously defined edges are reweighted.
%
%   A = GRASP_ADJACENCY_GAUSSIAN(graph, sigma) construct the adjacency
%       matrix A such that a_ij = exp(-d_ij² / (2 * sigma²)). Set sigma to
%       'auto' to use 1/3 of the average distance accross edges of the
%       graph (heuristic [2]). If the graph already has an adjacency
%       matrix, only distances corresponding to edges with non-zero weights
%       are used to compute the mean.
%
%   GRASP_ADJACENCY_GAUSSIAN(..., lambda) normalize weights with the power
%       lambda of each degree prior to computing the Gaussian kernel, i.e.
%       compute A=D^-lambda * A * D^-lambda [1].
%
%   [1] Graph Laplacians and their Convergence on Random Neighborhood 
%   Graphs. M. Hein, J.-Y. Audibert, and U. von Luxburg, 2007.
%
%   [2] Semi-Supervised Learning. O. Chapelle, B. Schölkopf, and A. Zien,
%   MIT Press, 2006.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2017-2020)
% Copyright Benjamin Girault, École Nationale de la Statistique et de
% l'Analyse de l'Information, Bruz, FRANCE (2020-2021)
% 
% benjamin.girault@ens-lyon.fr
% benjamin.girault@usc.edu
% benjamin.girault@ensai.fr

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

function A = grasp_adjacency_gaussian(graph, sigma, lambda)
    if nargin == 2
        lambda = 0;
    end
    
    N = size(graph.distances, 1);
    
    has_adja = false;
    if size(graph.A, 1) == N && size(graph.A, 2) == N
        has_adja = true;
        prev_A = (graph.A ~= 0);
    end
    
    if ischar(sigma)
        if strcmp(sigma, 'auto')
            tmp = graph.distances;
            if has_adja
                tmp = tmp(prev_A);
            end
            sigma = mean(tmp(:)) / 3;
        else
            error('Unrecogized paramater ''sigma'': it should be either a postive value or ''auto''.');
        end
    end
    if sigma <= 0
        error('Unrecogized paramater ''sigma'': it should be either a postive value or ''auto''.');
    end
    
    graph.A = exp(-graph.distances .^ 2 / (2 * sigma ^ 2)) - eye(N);
    if lambda ~= 0
        graph = grasp_adjacency_degreenorm(graph, lambda);
    end
    A = graph.A;
    if has_adja
        A = A .* prev_A;
    end
end
