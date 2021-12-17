%Shows the complete GFT using pseudo-time series.
%
%   GRASP_SHOW_TRANSFORM(ax_handle, graph) shows the GFT using the random
%       walk Laplacian 2nd eigenvector to order the nodes and embedded
%       them in 1D.
%
%   [embedding, clusters] = GRASP_SHOW_TRANSFORM(...) returns the 1D
%       embeddingt and the clusters used.
%
%   GRASP_SHOW_TRANSFORM(..., options) optional parameters:
%
%   options.transform_matrix: use the provided transform matrix instead of 
%       the GFT (default: [], for the GFT matrix). The transform matrix 
%       should be of size MxN, with N the number of nodes.
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
%       of node i (default: 1).
%
%   options.embedding: embedding of nodes in 1D (default: 0). Set to 0 
%       to use the random walk Laplacian embedding (default) or 1 for an
%       even embedding ordered according to the random walk Laplacian
%       embedding.
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
%   options.epsilon_support: threshold to consider a node in the support
%       of a mode or not (default: 0.05). Computed after scaling using
%       amplitude_normalization.
%
%   options.support_scatter_size: size of the dots for nodes in the 
%       support of the mode (default: 36). If dots are of variable width
%       (see options.support_scatter_mode), then this size corresponds to
%       an amplitude of 1.
%
%   options.support_scatter_mode: controls how the scatter plot of the
%       support is showing with respect to the amplitude and sign of the
%       mode. Default is constant size and black (default: 'constant').
%       Other possible values are 'var_size' for a variable width of the
%       dots, 'var_gray' for variable shades of gray, 'var_color' for 
%       variable color intensity (blue for negative, red for positive, 
%       white for zero), and 'var_width_color' for a combination of
%       'var_color' and 'var_width', or 'var_width_gray' for a combination 
%       of 'var_gray' and 'var_width'.
%
%   options.graph_signal_y_scheme: either place regularly space all the
%       graph signals of the transform ('regular', default), or use the
%       provided graph frequencies ('freq').
%
%   options.right_scale_size: how large is the right scale compared to the
%       left. 1 is equal size (default: 0.075).
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
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018-2019)
% Copyright Benjamin Girault, École Nationale de la Statistique et de
% l'Analyse de l'Information, Bruz, FRANCE (2020-2021)
% 
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

function [embedding, clusters] = grasp_show_transform(fig_handle, graph, varargin)
    %% Parameters
    default_param = struct(...
        'transform_matrix', [],...
        'highlight_entries', [],...
        'graph_frequencies', [],...
        'clusters', 1,...
        'embedding', 0,...
        'amplitude_scale', 1.5,...
        'amplitude_normalization', 'max_abs',...
        'epsilon_support', 0.05,...
        'support_scatter_size', 36,...
        'support_scatter_mode', 'constant',...
        'graph_signal_y_scheme', 'regular',...
        'right_scale_size', 0.075,...
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
        fprintf('Graph signal transform vizualization\n');
        fprintf('\tDOI: <to appear>\n');
    end
    
    %% Initializations
    clf(fig_handle);
    N = grasp_nb_nodes(graph);
    
    %% Transform matrix and graph frequencies
    if numel(options.transform_matrix) == 0 && size(graph.Finv, 1) ~= N
        error('Please compute the GFT prior to using this function (using grasp_eigendecomposition) or provide a transform using options.transform_matrix!');
    end
    
    transform_is_gft = false;
    if numel(options.transform_matrix) > 0
        % Provided matrix
        if size(options.transform_matrix, 2) ~= N
            error('options.transform_matrix should be of size MxN where N is the number of nodes!');
        end
        if size(options.transform_matrix, 1) ~= N
            if numel(options.graph_frequencies) ~= size(options.transform_matrix, 1)
                error('Please provide the correct number of graph frequencies with options.graph_frequencies!');
            end
            % frequencies ok
        else
            if numel(options.graph_frequencies) == 0 && numel(graph.eigvals) == N
                disp('Using graph.eigvals as graph frequencies.');
                options.graph_frequencies = graph.eigvals;
            elseif numel(options.graph_frequencies) ~= size(options.transform_matrix, 1)
                error('Please provide as many graph frequencies in options.graph_frequencies as lines in options.transform_matrix!');
            end
            % frequencies ok
        end
        if sum(size(graph.Q) == [N N]) == 2
            modes = graph.Q ^ (-1) * options.transform_matrix';
        else
            modes = options.transform_matrix';
        end
    else
        % GFT
        transform_is_gft = true;
        options.transform_matrix = graph.F;
        options.graph_frequencies = graph.eigvals;
        modes = graph.Finv;
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
    
    %% Clusters
    if numel(options.clusters) ~= N
        if numel(options.clusters) ~= 1
            error('''clusters'' options should be either the number of clusters to compute or a vector giving the cluster id of each node!');
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
    end
    if size(options.clusters, 2) > 1
        options.clusters = options.clusters';
    end
    
    %% Number of clusters
    num_clusters = numel(unique(options.clusters));
    
    %% Vertex embedding
    if numel(options.embedding) == 1
        disp('Computing node embedding using the second eigenvector of the Random Walk Laplacian...');
        warning('off', 'MATLAB:eigs:IllConditionedA');
        [V, ~] = eigs(grasp_laplacian_standard(graph), grasp_degrees(graph), 2, 0);
        warning('on', 'MATLAB:eigs:IllConditionedA');
        [~, options.ordering] = sort(V(:, 2));
        if options.embedding == 0
            options.embedding = V(:, 2);
        
            % Checking ordering is consistent with clusters
            x_coords = options.clusters(options.ordering);
            tmp2 = sort(x_coords(x_coords(2:end) - x_coords(1:(end - 1)) ~= 0));
            if sum(tmp2(2:end) - tmp2(1:(end - 1)) == 0) > 0
                error('Vertex embedding and clusters are inconsistent (split cluster(s) in the embedding).');
            end
        elseif options.embedding == 1
            disp('Computing regular node embedding...');
            data_points = ((1:N)' - 1) / (N - 1);
            if num_clusters == 1
                options.embedding(options.ordering) = data_points;
            else
                % In that case, we:
                % - Order according to clusters, and then the already computed ordering
                node_ordering_wrt_clusters = zeros(N, 1);
                % - Compute the induced regular embedding

                prev_node = 0;
                orig_ordering_inv(options.ordering) = (1:N)';
                for c = 1:num_clusters
                    cluster_mask = options.clusters == c;

                    orig_ordering = orig_ordering_inv(cluster_mask)';

                    % New ordering
                    remap = zeros(N, 1);
                    remap(sort(orig_ordering)) = (1:numel(orig_ordering))';
                    node_ordering_wrt_clusters(cluster_mask) = prev_node + remap(orig_ordering);

                    prev_node = prev_node + numel(orig_ordering);
                end

                % Ordering as required (instead of permutation)
                tmp(node_ordering_wrt_clusters) = (1:N)'; % inverse permutation
                options.embedding(tmp) = data_points;
                options.ordering = tmp;
            end
        else
            error('GraSP:ShowTransform:EmbeddingCode', 'Unknown embedding scheme.');
        end
    elseif numel(options.embedding) == N
        if size(options.embedding, 2) > 1
            options.embedding = options.embedding';
        end
        [~, options.ordering] = sort(options.embedding);
        
        % Checking ordering is consistent with clusters
        x_coords = options.clusters(options.ordering);
        tmp2 = sort(x_coords(x_coords(2:end) - x_coords(1:(end - 1)) ~= 0));
        if sum(tmp2(2:end) - tmp2(1:(end - 1)) == 0) > 0
            error('Vertex embedding and clusters are inconsistent (split cluster(s) in the embedding).');
        end
    else
        error('GraSP:ShowTransform:EmbeddingSize', 'Incorrect embedding size.');
    end
    embedding = options.embedding;
    
    %% Cluster ordering: relabel the clusters according to options.embedding
    cluster_mean = zeros(max(options.clusters), 1);
    for c = unique(options.clusters)'
        cluster_mean(c) = mean(options.embedding(options.clusters == c));
    end
    [~, IX] = sort(cluster_mean);
    tmp(IX) = (1:numel(cluster_mean));
    options.clusters = tmp(options.clusters)';
    clusters = options.clusters;
    
    %% Vertical lines to separate the clusters
    cluster_boundaries = [];
    if num_clusters > 1
        cluster_boundaries = zeros(num_clusters + 1, 1);
        cluster_boundaries(1) = min(options.embedding);
        cluster_boundaries(num_clusters + 1) = max(options.embedding);
        for c = 2:num_clusters
            cluster_boundaries(c) = (min(options.embedding(options.clusters == c)) + max(options.embedding(options.clusters == c - 1))) / 2;
        end
    end
    
    %% Frequencies y axis coordinates
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
    if numel(options.bands) > 0 && (size(options.bands_colors, 1) ~= size(options.bands, 1) || size(options.bands_colors, 2) < 3 || size(options.bands_colors, 2) > 4)
        warning('options.bands_colors is of an incorrect size!');
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
    cur_xlim = cur_xlim + (cur_xlim(2) - cur_xlim(1)) / 100 * [-1 1];
    plot(repmat(cur_xlim', 1, M),...
         repmat(spectral_spacing, 2, 1),...
         'Color', (1 - min_gray) * [1 1 1]);
    hold('on');
    final_yticks = get(gca, 'YTick');
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
    
    % Vertical bars to identify nodes
    yleftlimits = [min(spectral_spacing) max(spectral_spacing)] + cur_amplitude_scale * [-1 1];
    plot(repmat(options.embedding(options.ordering)', 2, 1), repmat(yleftlimits', 1, numel(options.embedding(options.ordering)), 1), 'Color', (1 - min_gray) * [1 1 1]);
    
    % Dots where nodes are in the support
    support_mask = abs(all_ordered_modes(:)) > options.epsilon_support;
    support_x = all_x(support_mask);
    support_y = all_series(support_mask);
    support_size = options.support_scatter_size;
    support_color = 'k';
    
    var_width = 0;
    var_color = 0;
    var_gray = 0;
    switch options.support_scatter_mode
        case 'constant'
        support_x = all_x(abs(all_ordered_modes(:)) > options.epsilon_support);
        support_y = all_series(abs(all_ordered_modes(:)) > options.epsilon_support);
        case 'var_width'
            var_width = 1;
        case 'var_gray'
            var_gray = 1;
        case 'var_color'
            var_color = 1;
        case 'var_width_gray'
            var_width = 1;
            var_gray = 1;
        case 'var_width_color'
            var_width = 1;
            var_color = 1;
        otherwise
            error('Unknown, or not implemented ''options.support_scatter_mode''');
    end
    if var_width
        support_size = options.support_scatter_size * abs(all_ordered_modes(support_mask));
    end
    if var_gray
        support_color = 1 - repmat(abs(all_ordered_modes(support_mask)), 1, 3);
    end
    if var_color
        mode_amplitude = abs(all_ordered_modes(support_mask));
        mode_sign = sign(all_ordered_modes(support_mask));
        support_color = ones(numel(mode_sign), 3);
        support_color(mode_sign > 0, 2) = 1 - mode_amplitude(mode_sign > 0);
        support_color(mode_sign > 0, 3) = support_color(mode_sign > 0, 2);
        support_color(mode_sign < 0, 1) = 1 - mode_amplitude(mode_sign < 0);
        support_color(mode_sign < 0, 2) = support_color(mode_sign < 0, 1);
    end
    if var_color || var_gray || var_width
        scatter(support_x,...
                support_y,...
                support_size,...
                support_color,...
                'filled',...
                'MarkerFaceColor', 'flat',...
                'MarkerEdgeColor', 'k');
    else
        scatter(support_x,...
                support_y,...
                support_size,...
                support_color,...
                '.');
    end
    
    % Highlighting nodes
    ordered_highlight = options.highlight_entries(:, options.ordering)';
    scatter(all_x(ordered_highlight > 0), all_series(ordered_highlight > 0), options.support_scatter_size, 'r');
    
    ylim(yleftlimits);
    set(gca, 'YTick', final_yticks);
        
    % Cluster Boundaries
    ylim_tmp = ylim();
    plot([cluster_boundaries(2:(end - 1)) cluster_boundaries(2:(end - 1))]', repmat(ylim_tmp, numel(cluster_boundaries) - 2, 1)', 'k', 'LineWidth', 1.5);
    
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
        alpha = 0.1;
        if size(options.bands_colors, 2) == 4
            alpha = options.bands_colors(l, 4);
        end
        patch([cur_xlim(1) cur_xlim cur_xlim(2)], [cur_band fliplr(cur_band)], options.bands_colors(l, 1:3), 'FaceAlpha', alpha, 'LineStyle', 'none');
    end
    hold('off');

    % Left label 
    if transform_is_gft
        ylabel('Graph frequency index $l$', 'Interpreter', 'Latex');
    else
        ylabel('Transform mode index $l$', 'Interpreter', 'Latex');
    end
    xlim(cur_xlim);
    final_xticks = get(gca, 'XTick');
  
    % Right plot
    delta_cur_lim = cur_xlim(2) - cur_xlim(1);
    new_xlim = cur_xlim;
    new_xlim(2) = cur_xlim(2) + options.right_scale_size * delta_cur_lim;

    % Frequencies
    yyaxis right
    freq_xs = [cur_xlim(2) cur_xlim(2) + .75 * options.right_scale_size * delta_cur_lim new_xlim(2)];
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
        
        alpha = 0.1;
        if size(options.bands_colors, 2) == 4
            alpha = options.bands_colors(l, 4);
        end
        patch([freq_xs fliplr(freq_xs)], [low_freq fliplr(high_freq)], options.bands_colors(l, 1:3), 'FaceAlpha', alpha, 'LineStyle', 'none');
    end
    hold('off');
    
    % Cleanup
    set(gca, 'YColor', zeros(1, 3));
    
%     freq_range = (options.graph_frequencies(M) - options.graph_frequencies(1));
%     vertical_clearance = options.amplitude_scale / (M - 1 + 2 * options.amplitude_scale);
%     vertical_clearance = vertical_clearance / (1 - 2 * vertical_clearance);
%     ylim(freq_range * [-vertical_clearance 1 + vertical_clearance])
    right_amplitude_scale = options.amplitude_scale * (options.graph_frequencies(M) - options.graph_frequencies(1)) / M;
    ylim(sort([options.graph_frequencies(1) options.graph_frequencies(M)]) + right_amplitude_scale * [-1 1]);
    
    ylabel('Graph frequency $\lambda_l$', 'Interpreter', 'Latex');
    
    xlim(new_xlim);
    set(gca, 'XTick', final_xticks);
    xlabel('Node embedding', 'Interpreter', 'Latex');
end
