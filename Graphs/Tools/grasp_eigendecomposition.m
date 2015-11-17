%Computes the eigendecomposition of a graph (Fourier operators).
%
%   graph = GRASP_EIGENDECOMPOSITION(graph) computes the Laplacian matrix,
%   and its eigendecomposition.
%
%   GRASP_EIGENDECOMPOSITION(..., options) optional parameters:
%
%   options.matrix: Use ne of the following matrices:
%       - 'std_lapl' (default): standard Laplacian eigendecomposition
%       - 'norm_lapl': normalized Laplacian eigendecomposition
%       - 'adja': adjacency matrix (graph shift)
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

function graph = grasp_eigendecomposition(graph, varargin)
    %% Parameters
    default_param = struct(...
        'matrix', 'std_lapl');
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    if isfield(options, 'laplacian')
        warning('DEPRECATED ''laplacian'' parameter! Use ''matrix'' instead.');
        switch options.laplacian
            case 'standard'
                options.matrix = 'std_lapl';
            case 'normalized'
                options.matrix = 'norm_lapl';
            otherwise
                error('Unrecognized Laplacian!');
        end
    end
    
    %% Checks
    if grasp_is_directed(graph)
        error('Error: directed graph not (yet) supported!');
    end
    
    %% Laplacian
    switch options.matrix
        case 'std_lapl'
            graph.L = grasp_laplacian_standard(graph);
            graph.fourier_version = 'standard laplacian';
            matrix = graph.L;
        case 'norm_lapl'
            graph.L = grasp_laplacian_normalized(graph);
            graph.fourier_version = 'normalized laplacian';
            matrix = graph.L;
        case 'adja'
            graph.L = graph.A;
            graph.fourier_version = 'graph shift';
            matrix = graph.L;
        otherwise
            error(['Unknown matrix parameter (' options.matrix ')!']);
    end
    
    %% Common part
    [graph.Finv, D] = eig(full(matrix));
    graph.eigvals = diag(D);
    
    % sorting...
    if strcmp(options.matrix, 'adja')
        [graph.eigvals, IX] = sort(graph.eigvals, 'descend');
    else
        [graph.eigvals, IX] = sort(graph.eigvals);
    end
    graph.Finv = graph.Finv(:, IX);
    graph.F = graph.Finv';  % graph.F is a unitary matrix
end