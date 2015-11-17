%Defines and returns an empty graph structure.
%
%   GRASP_STRUCT()
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

function g = grasp_struct
    g.A = 0;                      % Adjacency matrix
    g.A_layout = 0;               % Adjacency matrix to plot
    g.layout = 0;                 % Graph layout (in the 10x10 2D plane)
    g.distances = 0;              % Matrix of distances between any pair of nodes
    g.node_names = {};            % Cell of string for the names of nodes
    g.L = 0;                      % Laplacian matrix (std or norm)
    g.fourier_version = 'n.a.';   % Matrix used for the Fourier transform
    g.eigvals = 0;                % Eigenvalues of the graph
    g.F = 0;                      % Fourier matrix
    g.Finv = 0;                   % Invert Fourier matrix
    g.T = 0;                      % Translation operator
    g.background = '';            % Background file
end