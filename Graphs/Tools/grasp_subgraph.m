%Computes a subgraph.
%
%NOTE: this function only updates the standard fields of the graph
%structure (see GRASP_STRUCT).
%
%   g = GRASP_SUBGRAPH(g, vertex_set) remove the vertex not in vertex_set
%   from the graph g.
%
%   [..., mapping] = GRASP_SUBGRAPH(...) also outputs the mapping from the
%   old vertex index to the new vertex index (0 if not included).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2017-2018)
% 
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

function [g, mapping] = grasp_subgraph(g, vertex_set)
    Nold = grasp_nb_nodes(g);
    Nnew = numel(vertex_set);
    
    mapping = zeros(Nold, 1);
    mapping(vertex_set) = (1:Nnew)';
    
    if numel(g.A) > 1
        g.A = g.A(vertex_set, vertex_set);
    end
    if numel(g.A_layout) > 1
        g.A_layout = g.A_layout(vertex_set, vertex_set);
    end
    if numel(g.layout) > 1
        g.layout = g.layout(vertex_set, :);
    end
    if numel(g.distances) > 1
        g.distances = g.distances(vertex_set, vertex_set);
    end
    if numel(g.node_names) > 0
        g.node_names = g.node_names(vertex_set);
    end
    g.M = 0;
    g.Q = 0;
    g.Z = 0;
    g.fourier_version = 'n.a.';   % Matrix used for the Fourier transform
    g.eigvals = 0;                % Eigenvalues of the graph
    g.F = 0;                      % Fourier matrix
    g.Finv = 0;                   % Invert Fourier matrix
    g.T = 0;                      % Translation operator
end