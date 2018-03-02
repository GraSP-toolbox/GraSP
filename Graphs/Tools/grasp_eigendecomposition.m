%Computes the eigendecomposition of a graph (Fourier operators).
%
%   graph = GRASP_EIGENDECOMPOSITION(graph) computes the Laplacian matrix,
%   and its eigendecomposition.
%
%   GRASP_EIGENDECOMPOSITION(..., options) optional parameters:
%
%   options.matrix: Use one of the following matrices:
%       - 'std_lapl' (default): standard Laplacian eigendecomposition
%       - 'norm_lapl': normalized Laplacian eigendecomposition
%       - 'adja': adjacency matrix (graph shift)
%       - any valid NxN matlab matrix, with N the number of vertices
%       (matrix M in [1])
%   options.inner_product: Compute a familly of eigenvectors orthonormal
%       according to the given inner product (matrix Q in [1]).
%
%   [1] https://hal.archives-ouvertes.fr/hal-01708695
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Copyright Benjamin Girault, University of Southern California, USA
% (2017-2018).
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

function graph = grasp_eigendecomposition(graph, varargin)
    %% Parameters
    default_param = struct(...
        'matrix', 'std_lapl',...
        'inner_product', []);
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
        if ~strcmp(options.matrix, 'adja') || sum(size(options.matrix) ~= (grasp_nb_nodes(graph) * [1 1])) > 0
            error('Error: directed graph not (yet) supported!');
        else
            warning('The Fourier transform is most certainly not a unitary trasnform!');
        end
    end
    
    %% Laplacian
    if sum(size(options.matrix) ~= (grasp_nb_nodes(graph) * [1 1])) == 0
        options.fundamental_matrix = options.matrix;
        graph.L = options.matrix;
        options.matrix = 'custom';
    else
        if ~ischar(options.matrix)
            error('Unknown matrix parameter!');
        end
        switch options.matrix
            case 'std_lapl'
                graph.L = grasp_laplacian_standard(graph);
                graph.fourier_version = 'standard laplacian';
                options.fundamental_matrix = graph.L;
            case 'norm_lapl'
                graph.L = grasp_laplacian_normalized(graph);
                graph.fourier_version = 'normalized laplacian';
                options.fundamental_matrix = graph.L;
            case 'adja'
                graph.L = graph.A;
                graph.fourier_version = 'graph shift';
                options.fundamental_matrix = graph.L;
            otherwise
                error(['Unknown matrix parameter (' options.matrix ')!']);
        end
    end
    
    %% Inner Product
    if isempty(options.inner_product)
        graph.Q = speye(grasp_nb_nodes(graph));
    else
        graph.Q = options.inner_product;
    end
    
    %% Core
    if grasp_is_directed(graph)
        % Only the adjacency matrix is used here
        % We perform the Schur Block Diagonalization, more stable than the 
        % Jordan Normal Form [Girault, PhD Thesis, 2015]
        % The call to schur ensures a block trigonal matrix Tschur
        tmp = options.fundamental_matrix;
        if numel(options.inner_product) > 0
            tmp = (options.inner_product ^ -1) * tmp;
        end
        [V, T] = bdschur(tmp);
        [Y, Tblk] = schur(T, 'complex');
        graph.Finv = V * Y;
        graph.Tschur = Tblk;
        graph.eigvals = diag(graph.Tschur);
    else
        % Unitarily diagonalizable case
        if numel(options.inner_product) == 0
            [graph.Finv, D] = eig(full(options.fundamental_matrix));
        else
            if issparse(options.fundamental_matrix) || issparse(options.inner_product)
                % Matlab is stupid here: eig cannot be used to compute
                % eigenvalues and eigenvectors of sparse matrices (even in
                % the symmetric real case), yet eigs displays a warning in
                % that case that it falls back to eig, and make it full!
                %
                % We leave the warning on as this may be an issue when the
                % sparse matrix is large and its full version does not fit
                % into memory.
                opts.v0 = ones(grasp_nb_nodes(graph), 1);
%                 warning('off', 'MATLAB:eigs:TooManyRequestedEigsForRealSym');
                [graph.Finv, D] = eigs(options.fundamental_matrix, options.inner_product, grasp_nb_nodes(graph), 'sm', opts);
%                 warning('on', 'MATLAB:eigs:TooManyRequestedEigsForRealSym');
            else
                [graph.Finv, D] = eig(options.fundamental_matrix, options.inner_product);
            end
        end
        graph.eigvals = diag(D);
    end
    
    % Sorting...
    if strcmp(options.matrix, 'adja')
        if grasp_is_directed(graph)
            % Some rounding operation (we keep only 6 digits) to perform
            % sorting (otherwise equal eigvals may be swapped, and Tschur
            % does not remain upper triangular)
            [~, IX] = sort(round(abs(graph.eigvals), 5 + floor(log10(max(abs(graph.eigvals))))), 'descend');
            graph.eigvals = graph.eigvals(IX);
        else
            [graph.eigvals, IX] = sort(graph.eigvals, 'descend');
        end
    else
        [graph.eigvals, IX] = sort(graph.eigvals);
    end
    graph.Finv = graph.Finv(:, IX);
    if isfield(graph, 'Tschur')
        graph.Tschur = graph.Tschur(IX, IX);
    end
    if grasp_is_directed(graph)
        graph.F = graph.Finv ^ (-1);
    else
        % graph.Finv is a unitary matrix
        if numel(options.inner_product) == 0
            graph.F = graph.Finv';
        else
            graph.F = graph.Finv' * options.inner_product;
        end
    end
end