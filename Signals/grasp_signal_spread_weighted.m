%Computes the graph spread of a signal (Agaskar & Lu, IEEE Info. Th. 2015),
%but using the distance matrix provided (Pasdeloup, Alami, Gripon, Rabbat,
%EUSIPCO 2015).
%
%   T = GRASP_SIGNAL_SPREAD_WEIGHTED(graph) computes the weighted graph
%   spread of signal on graph, using graph.distances (graph.distances(i,j)
%   is the distance between nodes i and j).
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

function spread = grasp_signal_spread_weighted(graph, signal)
    %% Intialization
    if ~strcmp(graph.fourier_version, 'standard laplacian') && ~strcmp(graph.fourier_version, 'normalised laplacian')
        error('Error: The Fourier transform should be the result of the decomposition of a Laplacian matrix!');
    end
    
    %% Matrix of distances in the graph
    distances = graph.distances;
    
    %% Spread
    spread = 1 / norm(signal) ^ 2 * (distances .^ 2) * (abs(signal) .^ 2);
    spread = sqrt(min(spread));
end