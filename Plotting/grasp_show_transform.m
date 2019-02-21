%Shows the complete GFT using pseudo-time series.
%
%   GRASP_SHOW_TRANSFORM(ax_handle, graph) shows the GFT using the random
%       walk Laplacian 2nd eigenvector to order the vertices and embedded
%       them in 1D.
%
%   [embedding, clusters] = GRASP_SHOW_TRANSFORM(...) returns the 1D
%       embeddingt and the clusters used.
%
%   GRASP_SHOW_TRANSFORM(..., options) optional parameters:
%
%   options.transform_matrix: use the provided transform matrix instead of 
%       the GFT (default: [], for the GFT). The transform matrix should be
%       of size NxM, with N the number of vertices.
%
%   options.highlight_entries: 0/1 matrix of the same size than
%       options.transform_matrix indicating which entries to highlight in
%       the visualization (default: all 0).
%
%   options.graph_frequencies: assign a graph frequency to each of the M
%       atoms in the transform (default: graph.eigvals for the GFT).
%
%   options.clusters: either an integer for the number of clusters to
%       compute using spectral clustering based on the Random Walk 
%       Laplacian, or a vector such that clusters(i) is the cluster index 
%       of vertex i (default: 1).
%
%   options.ordering: ordering of vertices for the matrices (the vertex 
%       indices as they should appear in the embedding ordering(1) is the 
%       first vertex, and ordering(N) is the last). If ordering is
%       not of the right length, ordering will be computed from 
%       options.embedding if set or otherwise from sorting the second 
%       eigenvector of the Random Walk Laplacian (default: []).
%
%   options.embedding: embedding of vertices in 1D (default: [] for evenly
%       spread according to options.ordering). Set to 0 to use the random
%       walk Laplacian embedding.
%
%   options.amplitude_scale: scaling of the transform modes. Use the value 
%       1 for an optimal scaling without overlap, or 2 for an optimal 
%       scaling with exactly two modes overlapping (default: 1.5, i.e. some 
%       overlapping).
%
%   options.amplitude_normalization: 'l2', 'max_abs' or 'overall_max_abs'.
%       Use 'l2' to scale every transform mode to unit l2-norm, 'max_abs' 
%       to scale every transform mode to have unit maximum amplitude, and 
%       'overall_max_abs' to have the transform matrix maximum amplitude
%       equal to 1 (default: 'max_abs').
%
%   options.epsilon_support: threshold to consider a vertex in the support
%       of a mode or not (default: 0.05). Computed after scaling using
%       amplitude_normalization.
%
%   options.support_scatter_size: size of the dots for vertices in the 
%       support of the mode (default: 36).
%
%   options.graph_signal_y_scheme: either place regularly space all the
%       graph signals of the transform ('regular', default), or use the
%       provided graph frequencies ('freq').
%
%   options.bands: set of spectral bands to plot as subsets of graph 
%       frequencies: bands(i, 1) is the first graph frequency, and 
%       bands(i, 2) is the last for band i (default: [], i.e. no band).
%       This assumes that bands are ordered increasingly.
%
%   options.bands_colors: color for the shaded region of each band 
%       (default: [], i.e. default matlab colors).
%
%   options.verbose: whether or not to show the paper DOI (default: true).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018-2019)
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

function [embedding, clusters] = grasp_show_transform(fig_handle, graph, varargin)
    %% Parameters
    default_param = struct(...
        'transform_matrix', [],...
        'highlight_entries', [],...
        'graph_frequencies', [],...
        'clusters', 1,...
        'ordering', [],...
        'embedding', [],...
        'amplitude_scale', 1.5,...
        'amplitude_normalization', 'max_abs',...
        'epsilon_support', 0.05,...
        'support_scatter_size', 36,...
        'graph_signal_y_scheme', 'regular',...
        'bands', [],...
        'bands_colors', [],...
        'verbose', true);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Reference
    if options.verbose
        fprintf('\tDOI: <to appear>\n');
    end
    
    %% Initializations
    clf(fig_handle);
    N = grasp_nb_nodes(graph);
    
    %% Transform matrix and graph frequencies
    if numel(options.transform_matrix) == 0 && size(graph.Finv, 1) ~= N
        error('Please compute the GFT prior to using this function (using grasp_eigendecomposition) or provide a transform using options.transform_matrix!');
    end
    
    if numel(options.transform_matrix) > 0
        % Provided matrix
        if size(options.transform_matrix, 1) ~= N
            error('options.transform_matrix should be of size NxM where N is the number of vertices!');
        end
        if size(options.transform_matrix, 2) ~= N
            if numel(options.graph_frequencies) ~= size(options.transform_matrix, 2)
                error('Please provide the graph frequencies with options.graph_frequencies!');
            end
            % frequencies ok
        else
            if numel(options.graph_frequencies) == 0
                disp('Using graph.eigvals as graph frequencies.');
                options.graph_frequencies = graph.eigvals;
            elseif numel(options.graph_frequencies) ~= size(options.transform_matrix, 2)
                error('Please provide as many graph frequencies in options.graph_frequencies as columns in options.transform_matrix!');
            end
            % frequencies ok
        end
    else
        % GFT
        options.transform_matrix = graph.Finv;
        options.graph_frequencies = graph.eigvals;
    end
    if size(options.graph_frequencies, 2) > 1
        options.graph_frequencies = options.graph_frequencies';
    end
    
    %% Highlight Entries Matrix
    if numel(options.highlight_entries) == 0
        options.highlight_entries = zeros(size(options.transform_matrix));
    end
    if sum(size(options.highlight_entries) ~= size(options.transform_matrix)) > 0
        error('options.highlight_entries should have the same size than options.transform_matrix.');
    end
    
    %% Vertex ordering and embedding
    V = [];
    if numel(options.embedding) == 1 && options.embedding == 0
        if numel(V) == 0
            warning('off', 'MATLAB:eigs:IllConditionedA');
            [V, ~] = eigs(grasp_laplacian_standard(graph), grasp_degrees(graph), 2, 0);
            warning('on', 'MATLAB:eigs:IllConditionedA');
        end
        options.embedding = V(:, 2);
    end
    
    if numel(options.ordering) ~= N
        if numel(options.embedding) == N
            [~, options.ordering] = sort(options.embedding);
        else
            if numel(V) == 0
                warning('off', 'MATLAB:eigs:IllConditionedA');
                [V, ~] = eigs(grasp_laplacian_standard(graph), grasp_degrees(graph), 2, 0);
                warning('on', 'MATLAB:eigs:IllConditionedA');
            end
            disp('Computing vertex ordering using the second eigenvector of the Random Walk Laplacian...');
            [~, options.ordering] = sort(V(:, 2));
        end
    elseif numel(options.embedding) == N
        [~, tmp] = sort(options.embedding);
        if sum(arrayfun(@(i,j) double(i ~= j), options.ordering, tmp)) > 0
            error('options.ordering and options.embedding are not consistent!');
        end
    end
    if size(options.ordering, 2) > 1
        options.ordering = options.ordering';
    end
    
    %% Clusters
    if numel(options.clusters) ~= N
        if numel(options.clusters) ~= 1
            error('''clusters'' options should be either the number of clusters to compute or a vector giving the cluster id of each vertex!');
        end
        if options.clusters == 1
            options.clusters = ones(N, 1);
        elseif options.clusters < 1
            error('''clusters'' options should be at least 1!');
        else
            % Random Walk Laplacian clustering
            disp('Computing clusters using the Random Walk Clustering method...');
            warning('off', 'MATLAB:eigs:IllConditionedA');
            [V, ~] = eigs(grasp_laplacian_standard(graph), grasp_degrees(graph), round(options.clusters) + 1, 0);
            warning('on', 'MATLAB:eigs:IllConditionedA');
            options.clusters = kmeans(V, options.clusters);
        end
    else
        x_coords = options.clusters(options.ordering);
        tmp2 = sort(x_coords(x_coords(2:end) - x_coords(1:(end - 1)) ~= 0));
        if sum(tmp2(2:end) - tmp2(1:(end - 1)) == 0) > 0
            error('Vertex ordering and clusters are inconsistent.');
        end
    end
    
    %% Cluster ordering: relabel the clusters according to options.ordering
    cluster_mean = zeros(max(options.clusters), 1);
    tmp(options.ordering) = (1:N)';
    num_clusters = numel(cluster_mean);
    for c = 1:num_clusters
        cluster_mean(c) = mean(tmp(options.clusters == c));
    end
    [~, IX] = sort(cluster_mean);
    options.clusters = IX(options.clusters);
    
    clusters = options.clusters;
    
    %% Vertical lines to separate the clusters (only if using the regular
    % vertex embedding)
    cluster_boundaries = [];
    
    %% Vertex embedding
    if numel(options.embedding) ~= N
        disp('Computing regular vertex embedding...');
        data_points = ((1:N)' - 1) / (N - 1);
        if num_clusters == 1
            options.embedding(options.ordering) = data_points;
        else
            % In that case, we:
            % - Order according to clusters, and then the already computed ordering
            vertex_ordering_wrt_clusters = zeros(N, 1);
            % - Compute the induced regular embedding
            % - Compute boundaries of each cluster in the embedding
            clusters_first_and_last = zeros(num_clusters, 2);
            
            prev_vertex = 0;
            orig_ordering_inv(options.ordering) = (1:N)';
            for c = 1:num_clusters
                cluster_mask = options.clusters == c;
                
                orig_ordering = orig_ordering_inv(cluster_mask)';
                
                % New ordering
                remap = zeros(N, 1);
                remap(sort(orig_ordering)) = (1:numel(orig_ordering))';
                vertex_ordering_wrt_clusters(cluster_mask) = prev_vertex + remap(orig_ordering);
                
                prev_vertex = prev_vertex + numel(orig_ordering);
                
                % Boundaries of the cluster
                [~, IX] = min(vertex_ordering_wrt_clusters + (N * (1 - cluster_mask)));
                clusters_first_and_last(c, 1) = IX;
                [~, IX] = max(vertex_ordering_wrt_clusters .* cluster_mask);
                clusters_first_and_last(c, 2) = IX;
            end
            
            % Ordering as required (instead of permutation)
            tmp(vertex_ordering_wrt_clusters) = (1:N)'; % inverse permutation
            options.embedding(tmp) = data_points;
            options.ordering = tmp;
            
            % Boundaries of clusters in the embedded space
            tmp = clusters_first_and_last';
            cluster_boundaries(1) = options.embedding(tmp(1));
            cluster_boundaries(2:num_clusters) = mean([options.embedding(tmp(3:2:(end - 1))), options.embedding(tmp(2:2:(end - 2)))]);
            cluster_boundaries(num_clusters + 1) = options.embedding(tmp(end));
        end
    end
    if size(options.embedding, 2) > 1
        options.embedding = options.embedding';
    end
    embedding = options.embedding;
    
    %% Frequencies y axis coordinates
    modes = options.transform_matrix;
    M = size(modes, 2);
    switch options.graph_signal_y_scheme
        case 'regular'
            spectral_spacing  = 0:(M - 1);
        case 'freq'
            spectral_spacing  = options.graph_frequencies';
        otherwise
            error('Unknown parameter options.graph_signal_y_scheme.');
    end
    
    %% Bands
    if numel(options.bands) > 0
        if size(options.bands, 2) ~= 2
            error('options.bands should have two columns defining the spectral bands!');
        end
        options.bands_indices = horzcat(...
            arrayfun(@(lmin) find(options.graph_frequencies >= lmin, 1, 'first'), options.bands(:, 1), 'UniformOutput', false),...
            arrayfun(@(lmax) find(options.graph_frequencies <= lmax, 1, 'last'),  options.bands(:, 2), 'UniformOutput', false));
    else
        options.bands_indices = [];
    end
    
    %% Bands Colors
    if numel(options.bands) > 0 && sum(size(options.bands_colors) == [numel(options.bands) 2]) < 2
        options.bands_colors = [];
        while size(options.bands, 1) > size(options.bands_colors, 1)
            options.bands_colors = [options.bands_colors ; get(gca,'colororder')];
        end
    end
    
    %% Modes scaling
    switch options.amplitude_normalization
        case 'l2'
            modes = modes * diag(1 ./ arrayfun(@(i) norm(modes(:, i)), 1:size(modes, 2)));
        case 'max_abs'
            modes = modes * diag(1 ./ max(abs(modes)));
        case 'overall_max_abs'
            modes = modes / max(abs(modes(:)));
            % Nothing to do...
        otherwise
            error('unrecognized options.amplitude_normalization!');
    end
    
    %% Plotting
    min_gray = 0.1;
    
    cur_xlim = [min(options.embedding) max(options.embedding)];
    plot(repmat(cur_xlim', 1, M),...
         repmat(spectral_spacing, 2, 1),...
         'Color', (1 - min_gray) * [1 1 1]);
    hold('on');
    final_yticks = get(gca, 'YTick');Weights'
    final_yticks = final_yticks(final_yticks - floor(final_yticks) < 0.01);
    
    cur_amplitude_scale = options.amplitude_scale * (max(spectral_spacing) - min(spectral_spacing)) / M;
    cur_curve_scale = cur_amplitude_scale / 2;
    
    all_x             = repmat(options.embedding(options.ordering), 1, M);
    all_ordered_modes = modes(options.ordering, :);
    all_series        = cur_curve_scale * all_ordered_modes + repmat(spectral_spacing, N, 1);
    for l = 1:M
        y = all_ordered_modes(:, l);
        % Shades of a colormap version
%         colorindex = uint8((.2 + abs(y) / max(abs(y))) / 1.2 * 99) + 1; % TODO: change uint8 to some int
%         cdata = [uint8(mymap(101 - colorindex, :) * 255) uint8(ones(size(mymap(101 - colorindex, :), 1), 1))].';
        % Transparency shades version
%         coloralpha = uint8((min_gray + abs(y) / max(abs(y))) / (1 + min_gray) * 255);
        coloralpha = uint8((min_gray + abs(y)) / (1 + min_gray) * 255);
        cdata = [uint8(zeros(numel(coloralpha), 3)) coloralpha].';
        
        % Series plots (drawnow necessary to correctly shade the curves)
        ph = plot(all_x(:, l), all_series(:, l), '-k');
        drawnow;
        set(ph.Edge, 'ColorType', 'truecoloralpha', 'ColorBinding', 'interpolated', 'ColorData', cdata);
    end
    
    % Vertical bars to identify vertices
    yleftlimits = [min(spectral_spacing) max(spectral_spacing)] + cur_amplitude_scale * [-1 1];
    plot(repmat(options.embedding(options.ordering)', 2, 1), repmat(yleftlimits', 1, numel(options.embedding(options.ordering)), 1), 'Color', (1 - min_gray) * [1 1 1]);
    
    % Dots where vertices are in the support
    scatter(all_x(abs(all_ordered_modes(:)) > options.epsilon_support), all_series(abs(all_ordered_modes(:)) > options.epsilon_support), options.support_scatter_size, '.k');
    
    % Highlighting nodes
    ordered_highlight = options.highlight_entries(options.ordering, :);
    scatter(all_x(ordered_highlight > 0), all_series(ordered_highlight > 0), options.support_scatter_size, 'r');
    
    ylim(yleftlimits);
    set(gca, 'YTick', final_yticks);
        
    % Cluster Boundaries
    ylim_tmp = ylim();
    plot([cluster_boundaries(2:(end - 1)) cluster_boundaries(2:(end - 1))], repmat(ylim_tmp, numel(cluster_boundaries) - 2, 1), 'k', 'LineWidth', 1.5);
    
    % Bands
%     left_bands = zeros(size(options.bands_indices, 1), 2);
    %TODO: What if there is some overlap?
    left_bands = zeros(0, 2);
    for l = 1:size(options.bands_indices, 1)
        cur_band = options.bands_indices(l, :);
        if numel(cur_band{1}) == 0 || numel(cur_band{2}) == 0
            continue;
        end
        left_bands(l, :) = spectral_spacing(cell2mat(cur_band));
    end
    delta_left_band = (spectral_spacing(2) - spectral_spacing(1)) / 2;
    for l = 1:size(left_bands, 1)
        cur_band = left_bands(l, :);
        if cur_band(2) == 0 && cur_band(1) == 0
            continue;
        end
        if l > 1
            cur_band(1) = cur_band(1) - (cur_band(1) - left_bands(l - 1, 2)) / 2;
        else
            cur_band(1) = cur_band(1) - delta_left_band;
        end
        if l < size(left_bands, 1)
            cur_band(2) = cur_band(2) + (left_bands(l + 1, 1) - cur_band(2)) / 2;
        else
            cur_band(2) = cur_band(2) + delta_left_band;
        end
        patch([cur_xlim(1) cur_xlim cur_xlim(2)], [cur_band fliplr(cur_band)], options.bands_colors(l, :), 'FaceAlpha', 0.1, 'LineStyle', 'none');
    end
    hold('off');

    % Left label 
    ylabel('Graph Frequency index $l$', 'Interpreter', 'Latex');
    xlim(cur_xlim);
    final_xticks = get(gca, 'XTick');
  
    % Right plot
    delta_cur_lim = cur_xlim(2) - cur_xlim(1);
    new_xlim = cur_xlim;
    new_xlim(2) = cur_xlim(2) + .1 * delta_cur_lim;

    % Frequencies
    yyaxis right
    freq_xs = [cur_xlim(2) cur_xlim(2) + .075 * delta_cur_lim new_xlim(2)];
    freq_ys = [options.graph_frequencies(1) + (options.graph_frequencies(M) - options.graph_frequencies(1)) / (max(spectral_spacing) - min(spectral_spacing)) * spectral_spacing'...
               options.graph_frequencies...
               options.graph_frequencies];
    plot(repmat(freq_xs, M, 1)',...
         freq_ys',...
         '-', 'Color', [0.7 0.7 0.7])
    
    % Bands
    hold('on');
    for l = 1:size(options.bands_indices, 1)
        cur_band = options.bands_indices(l, :);
        if numel(cur_band{1}) == 0 || numel(cur_band{2}) == 0
            continue;
        end
        cur_band = cell2mat(cur_band);
        
        if cur_band(1) > 1
            low_freq = [mean([freq_ys(cur_band(1), 1) freq_ys(cur_band(1) - 1, 1)]) options.bands(l, 1) options.bands(l, 1)];
        else
            low_freq = (freq_ys(cur_band(1), 1) - freq_ys(2, 1) / 2) * ones(1, 3);
        end
        
        if cur_band(2) < M
            high_freq = [mean([freq_ys(cur_band(2), 1) freq_ys(cur_band(2) + 1, 1)]) options.bands(l, 2) options.bands(l, 2)];
        else
            high_freq = (freq_ys(cur_band(2), 1) + freq_ys(2, 1) / 2) * ones(1, 3);
        end
        
        patch([freq_xs fliplr(freq_xs)], [low_freq fliplr(high_freq)], options.bands_colors(l, :), 'FaceAlpha', 0.1, 'LineStyle', 'none');
    end
    hold('off');
    
    % Cleanup
    set(gca, 'YColor', [0.7 0.7 0.7]);
    
%     freq_range = (options.graph_frequencies(M) - options.graph_frequencies(1));
%     vertical_clearance = options.amplitude_scale / (M - 1 + 2 * options.amplitude_scale);
%     vertical_clearance = vertical_clearance / (1 - 2 * vertical_clearance);
%     ylim(freq_range * [-vertical_clearance 1 + vertical_clearance])
    right_amplitude_scale = options.amplitude_scale * (options.graph_frequencies(M) - options.graph_frequencies(1)) / M;
    ylim([options.graph_frequencies(1) options.graph_frequencies(M)] + right_amplitude_scale * [-1 1]);
    
    ylabel('Graph Frequency $\lambda_l$', 'Interpreter', 'Latex');
    
    xlim(new_xlim);
    set(gca, 'XTick', final_xticks);
    xlabel('Vertex embedding', 'Interpreter', 'Latex');
end
