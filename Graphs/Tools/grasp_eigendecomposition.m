%Computes the eigendecomposition of a graph (Fourier operators).
%
%   graph = GRASP_EIGENDECOMPOSITION(graph) computes the Laplacian matrix,
%       and its eigendecomposition.
%
%   GRASP_EIGENDECOMPOSITION(..., options) optional parameters:
%
%   options.matrix: Use one of the following matrices:
%       - 'std_lapl' (default): standard Laplacian eigendecomposition
%       - 'norm_lapl': normalized Laplacian eigendecomposition
%       - 'adja': adjacency matrix (graph shift)
%       - any valid NxN matlab matrix, with N the number of vertices
%       (matrix M in [1])
%
%   options.inner_product: Compute a familly of eigenvectors orthonormal
%       according to the given inner product (matrix Q in [1]), or one of
%       the named inner products:
%       - 'dotprod': dot product (default)
%       - 'degree': degree matrix of graph
%
%   [1] https://hal.archives-ouvertes.fr/hal-01708695
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2017-2019)
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

function graph = grasp_eigendecomposition(graph, varargin)
    %% Parameters
    default_param = struct(...
        'matrix', 'std_lapl',...
        'inner_product', [],...
        'verbose', true);
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    if isfield(options, 'laplacian')
        warning('GraSP:Deprecated', 'DEPRECATED ''laplacian'' parameter! Use ''matrix'' instead.');
        switch options.laplacian
            case 'standard'
                options.matrix = 'std_lapl';
            case 'normalized'
                options.matrix = 'norm_lapl';
            otherwise
                error('GraSP:Eigendecomposition:UnknownLaplacian', 'Unrecognized Laplacian!');
        end
    end
    
    %% Checks
    is_invalid_matrix = ~ischar(options.matrix) && sum(size(options.matrix) ~= (grasp_nb_nodes(graph) * [1 1])) > 0;
    if is_invalid_matrix
        error('GraSP:Eigendecomposition:MatrixIncorrectSize', 'Error: options.matrix is of incorrect size!');
    end
    
    if grasp_is_directed(graph)
        is_adja = ischar(options.matrix) && strcmp(options.matrix, 'adja');
        is_nonhermitian = ~ischar(options.matrix) && ~ishermitian(options.matrix);
        if is_adja
            warning('GraSP:Eigendecomposition:GFTNonUnitary', 'The Fourier transform is most certainly not a unitary transform!');
        elseif is_nonhermitian
            error('GraSP:Eigendecomposition:MatrixNonHermitian', 'Error: Non Hermitian matrix which is not the adjacency matrix is not supported!');
        end
    end
    
    %% Quadratic Graph Variation Operator
    if sum(size(options.matrix) ~= (grasp_nb_nodes(graph) * [1 1])) == 0
        options.variation_matrix = options.matrix;
        graph.M = options.matrix;
        options.matrix = 'custom';
        if ~ismatrix(graph.M) || size(graph.M, 1) ~= size(graph.M, 2) || size(graph.M, 1) ~= grasp_nb_nodes(graph)
            error('GraSP:Eigendecomposition:MatrixNonSquare', 'options.matrix should be a square matrix of the same size than the graph!');
        end
        if ishermitian(graph.M)
            [~, p] = chol(graph.M + 1e-15 * eye(grasp_nb_nodes(graph))); % Trick to quickly test whether M is semi-definite positive (up to an 1e-15 error)
            if p == 0
                graph.fourier_version = 'irregularity-aware';
            else
                warning('GraSP:Eigendecomposition:GFTDegenerateVariation', 'Degenerate irregularity-aware GFT: option.matrix is not semi-definite positive!');
                graph.fourier_version = 'degenerate irregularity-aware';
            end
        else
            % Non Hermitian case => we don't have a quadratic graph variation.
            graph.fourier_version = 'diagonalization';
        end
    elseif ischar(options.matrix)
        dotprod_cond = isempty(options.inner_product);
        dotprod_cond = dotprod_cond || (ischar(options.inner_product) && strcmp(options.inner_product, 'dotprod'));
        dotprod_cond = dotprod_cond || (isnumeric(options.inner_product) && nnz(options.inner_product - speye(grasp_nb_nodes(graph))) == numel(graph.A));
        switch options.matrix
            case 'std_lapl'
                graph.M = grasp_laplacian_standard(graph);
                if dotprod_cond
                    graph.fourier_version = 'standard laplacian';
                elseif ischar(options.inner_product) && strcmp(options.inner_product, 'degree')
                    graph.fourier_version = 'random walk laplacian';
                else
                    graph.fourier_version = 'irregularity-aware';
                end
                options.variation_matrix = graph.M;
            case 'norm_lapl'
                graph.M = grasp_laplacian_normalized(graph);
                if dotprod_cond
                    graph.fourier_version = 'normalized laplacian';
                else
                    graph.fourier_version = 'irregularity-aware';
                end
                options.variation_matrix = graph.M;
            case 'adja'
                graph.M = graph.A;
                if dotprod_cond
                    graph.fourier_version = 'graph shift';
                else
                    graph.fourier_version = 'irregularity-aware';
                end
                options.variation_matrix = graph.M;
            otherwise
                error('GraSP:Eigendecomposition:MatrixUnknown', ['Unknown matrix parameter (' options.matrix ')!']);
        end
    else
        error('GraSP:Eigendecomposition:MatrixUnknown', 'Unknown matrix parameter!');
    end
    
    %% Graph Signals Inner Product
    if isempty(options.inner_product)
        graph.Q = speye(grasp_nb_nodes(graph));
    elseif isnumeric(options.inner_product)
        if ~ishermitian(options.inner_product)
            error('GraSP:Eigendecomposition:NonHermitianInnerProduct', 'options.inner_product needs to be Hermitian!');
        end
        [~, p] = chol(options.inner_product); % Trick to quickly test whether M is definite positive
        if p ~= 0
            error('GraSP:Eigendecomposition:NonPositiveInnerProduct', 'options.inner_product needs to be a positive matrix!');
        end
        graph.Q = options.inner_product;
    elseif ischar(options.inner_product)
        switch options.inner_product
            case 'dotprod'
                graph.Q = speye(grasp_nb_nodes(graph));
                options.inner_product = [];
            case 'degree'
                if grasp_is_directed(graph)
                    if strcmp(options.matrix, 'std_lapl')
                        graph.Q = diag(diag(graph.M));
                    else
                        error('GraSP:Eigendecomposition:InnerProductDegreeMatrixStdLaplWithDirectedGraph', 'options.inner_product set to ''degree'' without using the standard Laplacian and with a directed graph (unsupported)!');
                    end
                else
                    graph.Q = grasp_degrees(graph);
                end
                options.inner_product = graph.Q;
            otherwise
                error('GraSP:Eigendecomposition:InnerProductUnknown', ['Unknown inner product type ''' options.inner_product '''!']);
        end
    else
        error('GraSP:Eigendecomposition:InnerProductUnknown', 'options.inner_product should be a matrix or string!');
    end
    
    %% Core
    if grasp_is_directed(graph) && strcmp(options.matrix, 'adja')
        % The adjacency matrix is used here
        % We perform the Schur Block Diagonalization, more stable than the 
        % Jordan Normal Form [Girault, PhD Thesis, 2015]
        % The call to schur ensures a block trigonal matrix Tschur
        if options.verbose
            fprintf('Graph Shift GFT\n')
            fprintf('\tDOI: https://doi.org/10.1109/TSP.2013.2238935\n');
            fprintf('Using Schur block diagonalization instead of Jordan decomposition\n');
            fprintf('\tHAL: https://tel.archives-ouvertes.fr/tel-01256044\n');
        end

        tmp = options.variation_matrix;
        if numel(options.inner_product) > 0
            tmp = (options.inner_product ^ -1) * tmp;
        end
        [V, T] = bdschur(tmp);
        [Y, Tblk] = schur(T, 'complex');
        graph.Finv = V * Y;
        graph.Tschur = Tblk;
        graph.eigvals = diag(graph.Tschur);
        graph.Z = options.variation_matrix;
    else
        % Unitarily diagonalizable case
        if numel(options.inner_product) == 0
            if options.verbose
                fprintf('Standard diagonalization\n');
                fprintf('\tDOI: https://doi.org/10.1109/MSP.2012.2235192\n');
            end
            
            [graph.Finv, D] = eig(full(options.variation_matrix));
            graph.Z = options.variation_matrix;
        else
            if options.verbose
                fprintf('Irregularity-Aware GFT\n');
                fprintf('\tDOI: https://doi.org/10.1109/TSP.2018.2870386\n');
            end
            
            if issparse(options.variation_matrix) || issparse(options.inner_product)
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
                [graph.Finv, D] = eigs(options.variation_matrix, options.inner_product, grasp_nb_nodes(graph), 'sm', opts);
%                 warning('on', 'MATLAB:eigs:TooManyRequestedEigsForRealSym');
            else
                [graph.Finv, D] = eig(options.variation_matrix, options.inner_product);
            end
            graph.Z = options.inner_product ^ (-1) * options.variation_matrix;
        end
        graph.eigvals = diag(D);
    end
    
    %% Sorting of graph frequencies...
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
        if isempty(options.inner_product)
            graph.F = graph.Finv';
        else
            graph.F = graph.Finv' * options.inner_product;
        end
    end
end
