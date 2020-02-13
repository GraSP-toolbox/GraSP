%Plots a graph and a graph signal. See the options below for details on all
%the capabilities.
%
% Arrow definition (directed edges)
%
%An arrow is composed of a tail and two heads, one on each side. The 
%length of each of these heads is a percentage head_proportion of the
%length of the tail.
%The ends of the heads are such that they are always at the same distance
%of the tail, and form an angle of 45° with the tail for the shortest
%arrow.
%The tail of the arrow does not reach exactly the center of the node
%on its tip (for aestetic reasons), but is reduced by a quantity
%tail_back_quantity.
%
%   GRASP_SHOW_GRAPH(axis_handle, graph) plots the graph on the axis
%       pointed by axis_handle. Uses the field graph.show_graph_options for
%       default plotting values of GRASP_SHOW_GRAPH.
%
%   [nodes_handle, edges_handle] = GRASP_SHOW_GRAPH(...) returns the handle
%       to the nodes graphics objects (scatter), and edge graphics object. 
%       For the edges, it returns an array of 3 elements: the non directed 
%       edges (plot), the arrow heads (fill), and the arrow tails (plot)
%
%   GRASP_SHOW_GRAPH(..., options) optional parameters:
%
%               ***Plot***
%
%   options.background: plot over a background image (default: set by the
%       input graph).
%
%   options.background_boundaries: mapping between the boundaries of the
%       images and the layout coordinate system (default: same as the
%       layout_boundaries, see below).
%
%   options.background_grayscale: whether the background should be
%       converted to grayscale (default: true).
%
%   options.layout_boundaries: 2 rows "x", "y" (or 3 for 3D layout, "x", 
%       "y", "z") of with boundaries to plot the graph [min max] (default: 
%       [0 10] for each dimension). Set to a negative value for a automatic
%       computation of the boundaries with a 5% margin, or any positive
%       value in the interval [0,1) to set the margin size (0.05 = 5%).
%
%   options.axis_style: see AXIS style input property (default: equal).
%
%   options.viewpoint3D: arguments to give to VIEW to set the camera
%       position for a 3D plot.
%
%               ***Vertex Labels***
%
%   options.node_text: use the provided cell to plot also a label
%       associated to each node (disabled: cell(0) (default), string cell
%       with one element per node: custom labels, 'ID': id of the node,
%       anything else: graph.node_names if set, or id of the node
%       otherwise.)
%
%   options.node_text_shift: shift the node label relative to the node
%       (default: no shift).
%
%   options.node_text_fontsize: Font size of the node labels (default: 20).
%
%   options.node_text_background_color: the color of the background
%       underneath the node label (default: white, use 'none' for 
%       transparent).
%
%   options.node_text_background_edge: the color of the background edge
%       underneath the node label (default: black, use 'none' to disable).
%
%               ***Vertex Values***
%
%   options.node_values: plot the values on the nodes using colors
%       (default: no values).
%
%   options.value_scale: use the color scale provided (default: [min max]
%       of node_values).
%
%   options.color_map: use the provided named color map as defined in 
%       matlab COLORMAP function, or 'blue_red' for the one set by 
%       GRASP_SET_BLUE_RED_COLORMAP (default: 'default' or 'blue_red' if 
%       value_scale is [-x x] for some value x).
%
%   options.show_colorbar: show the color bar on the right (default:
%       false).
%
%   options.node_display_size: use the size provided for the nodes
%       (default: 200).
%
%   options.node_marker_edge_width: each node is represented by a marker,
%       this options sets the width of the edge of the marker, 0 disables
%       that edge (default: 1).
%
%               ***Highlighted Nodes***
%
%   options.highlight_nodes: draw a bigger circle around the given nodes.
%
%   options.highlight_nodes_size: size of the circle (default:
%       2 * node_display_size).
%
%   options.highlight_nodes_width: width of the circle line (default: 0.5).
%
%   options.highlight_nodes_color: color for the circle line (default:
%       red).
%
%               ***Edges***
%
%   options.show_edges: whether or not to show edges of the graph
%       (default: true if there is no background, false otherwise).
%
%   options.edge_color: use the provided color for the edges (default:
%       [0 0 0] ; use [0.9 0.9 0.9] for a light gray not visible on a
%       print).
%
%   options.edge_colormap: use a colormap to draw the edges instead of a
%       fixed color (edge_color becomes then ineffective). An empty string
%       disables this feature. Values associated to each edges are obtained
%       through the matrix graph.A_layout, such that the color may
%       correspond to something different than the edge weight.
%
%   options.edge_colorbar: whether or not to show a colorbar on the bottom
%       of the plot to display the color mapping to weight values for the
%       edges (default: false).
%
%   options.edge_color_scale: together with edge_colormap, set the
%       weight boundaries in the colormap.
%
%   options.edge_thickness: use the thickness provided to draw edges or
%       directed edges (default: 0.5).
%
%               ***Directed Edges***
%
%   options.arrow_max_tip_back_fraction: the tip of the arrow is pulled
%       slightly back. Fraction of the the edge length that the arrow tip
%       is actually pulled back (default: 5%).
%
%   options.arrow_max_head_fraction: maximum fraction of the edge length
%       dedicated to the arrow head (default: 70%). This is a limit for
%       edges of 0 length.
%
%   options.arrow_width_screen_fraction: width of the arrow head as a
%       fraction of the axes boundaries diagonal length (default: 0.5%).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016-2019)
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

function [nodes_handle, edges_handle] = grasp_show_graph(axis_handle, input_graph, varargin)
    %% Parameters
    default_param = struct(...
        'background', input_graph.background,...
        'background_boundaries', 0,...
        'background_grayscale', true,...
        'layout_boundaries', [],...
        'axis_style', 'equal',...
        'viewpoint3D', [20 45],...
        'node_values', zeros(grasp_nb_nodes(input_graph), 1),...
        'value_scale', 0,...
        'highlight_nodes', [],...
        'highlight_nodes_size', -1,...
        'highlight_nodes_width', 0.5,...
        'highlight_nodes_color', 'r',...
        'color_map', 'default',...
        'show_colorbar', false,...
        'node_display_size', 200,...
        'node_marker_edge_width', 1,...
        'node_text', [],...
        'node_text_shift', 0,...
        'node_text_fontsize', 20,...
        'node_text_background_color', 'white',...
        'node_text_background_edge', 'black',...
        'show_edges', true,...
        'edge_color', [0 0 0],...
        'edge_colormap', '',...
        'edge_colorbar', false,...
        'edge_color_scale', 0,...
        'edge_thickness', 0.5,...
        'arrow_max_tip_back_fraction', 0.02,...
        'arrow_max_head_fraction', 0.7,...
        'arrow_width_screen_fraction', 0.005);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    if isfield(input_graph, 'show_graph_options')
        default_param = grasp_merge_structs(default_param, input_graph.show_graph_options);
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Intialization
    N = grasp_nb_nodes(input_graph);
    A_layout = input_graph.A;
    if isfield(input_graph, 'A_layout') && numel(input_graph.A_layout) == N ^ 2
        A_layout = input_graph.A_layout;
        disp('Using adjacency matrix provided by ''A_layout''.');
    end
    
    %% Highlighted nodes default values
    if options.highlight_nodes_size == -1
        options.highlight_nodes_size = 2 * options.node_display_size;
    end

    %% Layout
    if ~isfield(input_graph, 'layout') || size(input_graph.layout, 1) ~= N
        warning('GraSP:MissingLayout', 'Computing a spectral layout. Consider defining it manually beforehand for efficiency.');
        input_graph.layout = grasp_layout_spectral(input_graph);
    end
    
    %% Setting the axes
    
    if isempty(options.layout_boundaries)
        options.layout_boundaries = kron([0 10], ones(size(input_graph.layout, 2), 1));
    elseif numel(options.layout_boundaries) == 1
        if options.layout_boundaries < 0
            margin = 0.05; % 5% margin
        else
            margin = options.layout_boundaries;
        end
        options.layout_boundaries = [min(input_graph.layout)' max(input_graph.layout)'];
        boundaries_size = options.layout_boundaries(:, 2) - options.layout_boundaries(:, 1);
        margin_size = margin * boundaries_size;
        if margin_size(1) == 0
            if margin_size(2) == 0
                margin_size = [1 ; 1];
            else
                margin_size(1) = margin_size(2);
            end
        end
        if margin_size(2) == 0
            margin_size(2) = margin_size(1);
        end
        
        options.layout_boundaries(:, 1) = options.layout_boundaries(:, 1) - margin_size;
        options.layout_boundaries(:, 2) = options.layout_boundaries(:, 2) + margin_size;
    elseif size(options.layout_boundaries, 2) ~= 2 || size(options.layout_boundaries, 1) < 2 || size(options.layout_boundaries, 1) > 3
        error('''layout_boundaries'' should be a 2x2 or a 3x2 matrix!');
    end
    
    prev_hold = get(axis_handle, 'NextPlot');
    if strcmp(prev_hold, 'replace')
        cla(axis_handle);
        % Remove any remaining edge colorbar
        for ch = findall(axis_handle.Parent.Children, 'type', 'Axes')'
            if ch == axis_handle
                continue;
            end
            if sum(abs(ch.Position - axis_handle.Position)) == 0 && strcmp(ch.Title.String, 'edge_colorbar')
                delete(ch);
            end
        end
    end
    hold(axis_handle, 'on');
    axis(options.axis_style);
    
    % Setting boundaries now for optimisation
    axis(axis_handle, [options.layout_boundaries(1, :) options.layout_boundaries(2, :)]);
    set(axis_handle, 'XTick', [], 'YTick', []);
    if size(input_graph.layout, 2) == 3  % 3D layout
        zlim(axis_handle, options.layout_boundaries(3, :));
        view(axis_handle, options.viewpoint3D);
        set(axis_handle, 'ZTick', []);
%         camzoom(axis_handle, 1.5);
    end
    set(axis_handle, 'Box', 'on');
    if ~isempty(options.edge_colormap)
        colormap(axis_handle, options.edge_colormap);
        options.edge_colormap = colormap(axis_handle);
        if numel(options.edge_color_scale) ~= 2
            options.edge_color_scale = [...
                full(min(A_layout(A_layout ~= 0)))...
                full(max(A_layout(A_layout ~= 0)))];
        end
    end
    
    if strcmp(options.color_map, 'default') && numel(options.value_scale) == 2 && options.value_scale(1) + options.value_scale(2) == 0
        grasp_set_blue_red_colormap(axis_handle);
    elseif strcmp(options.color_map, 'blue_red')
        grasp_set_blue_red_colormap(axis_handle);
    else
        colormap(axis_handle, options.color_map);
    end
    
    %% Background
    % Getting rid of any remaining background image
    if strcmp(prev_hold, 'replace')
        for ch = findall(axis_handle, 'Type', 'Image')'
            delete(ch);
        end
    end
    % Then adding the new (if required)
    if ~isempty(options.background)
        % Background file or image?
        if ischar(options.background)
            [img, map, alpha] = imread(options.background);
        else
            img = options.background;
            map = [];
            alpha = [];
        end
        % Converting RGB to data used by imagesc, possibly in grayscale
        if ~isempty(map)
            if options.background_grayscale
                map = rgb2gray(map);
            end
            img = ind2rgb(img, map);
        elseif options.background_grayscale
            img = rgb2gray(img);
            img(:, :, 2) = img(:, :, 1);
            img(:, :, 3) = img(:, :, 1);
        end
        % Flipping to get the origin on the bottom left
        for i = 1:size(img, 3)
            img(:, :, i) = flipud(img(:, :, i));
        end
        % Bugfix for GNU Octave ()
        if grasp_is_octave && numel(alpha) > 0
            warning('GraSP:OctaveBugfixTrick', 'Blending alpha channel!');
            alpha_mask = double(flipud(alpha)) / 255;
            for i = 1:size(img, 3)
                img(:, :, i) = uint8(alpha_mask .* img(:, :, i) + (1 - alpha_mask) * 255);
            end
            alpha = [];
        end
        % And finally, showing the image
        if numel(options.background_boundaries) == 4 && size(options.background_boundaries, 1) == 2
            imh = imagesc(options.background_boundaries(1, :), options.background_boundaries(2, :), img);
        else
            imh = imagesc(options.layout_boundaries(1, :), options.layout_boundaries(2, :), img);
        end
        if numel(alpha) > 0
            set(imh, 'AlphaData', flipud(alpha));
        end
    end
    
    %% Color scale
    if numel(options.value_scale) ~= 2 && numel(options.node_values) ~= 0
        options.value_scale = [min(options.node_values) max(options.node_values)];
    end

    if numel(options.value_scale) == 2
        if options.value_scale(1) > options.value_scale(2)
            options.value_scale = options.value_scale([2 1]);
        end
        if options.value_scale(1) == options.value_scale(2) || isnan(options.value_scale(1)) || isnan(options.value_scale(2))
            options.value_scale = [0 1];
        end
        caxis(axis_handle, options.value_scale);
    end
    
    %% Edges
    edges_handle = cell(3, 1);
    if options.show_edges
        % Select undirected edges
        adja = A_layout ~= 0;
        adja = adja - diag(diag(adja)); % remove the diagonal
        undir = adja + adja' > 1;
        [rows, cols] = find(undir); % non directed edges
        % Plot them
        if size(input_graph.layout, 2) == 2
            edges_handle{1} = plot(axis_handle,...
                                  [input_graph.layout(rows', 1)' ; input_graph.layout(cols', 1)'],...
                                  [input_graph.layout(rows', 2)' ; input_graph.layout(cols', 2)'],...
                                  '-',...
                                  'linewidth', options.edge_thickness,...
                                  'Color', options.edge_color);
        else
            set(gcf,'CurrentAxes',axis_handle)
            edges_handle{1} = plot3([input_graph.layout(rows', 1)' ; input_graph.layout(cols', 1)'],...
                                    [input_graph.layout(rows', 2)' ; input_graph.layout(cols', 2)'],...
                                    [input_graph.layout(rows', 3)' ; input_graph.layout(cols', 3)'],...
                                    '-',...
                                    'linewidth', options.edge_thickness,...
                                    'Color', options.edge_color);
        end
        
        % Adjust edge color if necessary
        if ~isempty(options.edge_colormap)
            flipped_cm = flipud(options.edge_colormap);
            w = A_layout(sub2ind(size(A_layout), rows, cols));
            min_w = options.edge_color_scale(1);
            span_w = options.edge_color_scale(2) - options.edge_color_scale(1);
            colors = (w - min_w) / (span_w);
            colors = 1 + (1 - min(max(colors, 0), 1)) * (size(flipped_cm, 1) - 1);
            colors = flipped_cm(floor(colors), :);
            for i = 1:numel(edges_handle{1})
                if verLessThan('matlab', '8.4.0')
                    set(edges_handle{1}(i), 'Color', colors(i, :));
                else
                    edges_handle{1}(i).Color = colors(i, :);
                end
            end
        end

        % Select directed edges
        [rows, cols] = find(adja - undir); % directed edges
        if numel(rows) > 0
            % Edge lengths
            x_orig = input_graph.layout(rows', 1)';
            y_orig = input_graph.layout(rows', 2)';
            x_targ = input_graph.layout(cols', 1)';
            y_targ = input_graph.layout(cols', 2)';
            dx = x_targ - x_orig;
            dy = y_targ - y_orig;
            
            % Screen length of the arrow tail (assuming square axes)
            boundaries_size = options.layout_boundaries(:, 2) - options.layout_boundaries(:, 1);
            edge_screen_size = sqrt((dx / boundaries_size(1)) .^ 2 + (dy / boundaries_size(2)) .^ 2);
            
            % Computing head and tip back fractions
            base_fraction = 1 ./ (1 + 30 * edge_screen_size);
            head_fraction = options.arrow_max_head_fraction * base_fraction;
            tip_back_fraction = options.arrow_max_tip_back_fraction * (1 - base_fraction);
            
            % Final arrow lengths
            x_targ = x_targ - tip_back_fraction .* dx;
            y_targ = y_targ - tip_back_fraction .* dy;
            
            % Direction of the arrow
            main_vector_direction = -[x_targ - x_orig ; y_targ - y_orig];
            main_vector_direction = cell2mat(arrayfun(@(i) main_vector_direction(:, i) / norm(main_vector_direction(:, i)), (1:size(main_vector_direction, 2)), 'UniformOutput', false));
            
            % Head computations
            arrow_dx     = head_fraction .* dx;
            arrow_dy     = head_fraction .* dy;
            arrow_length = sqrt(arrow_dx .^ 2 + arrow_dy .^ 2);
            arrow_width  = options.arrow_width_screen_fraction * norm(boundaries_size);
            
            arrow_angle = atan(arrow_width ./ arrow_length);
            arrow_side_length = sqrt((arrow_width / 2) .^ 2 + (arrow_length) .^ 2);

            % Computing the head (for filling a triangle) and the tail (for
            % the edge plotting)
            M = numel(x_orig);
            headsX = NaN * zeros(3, M);
            headsY = NaN * zeros(3, M);
            tailsX = NaN * zeros(1, 3 * M);
            tailsY = NaN * zeros(1, 3 * M);
            for i = 1:M
                rotation_1 = [cos(arrow_angle(i) / 2) -sin(arrow_angle(i) / 2) ; sin(arrow_angle(i) / 2) cos(arrow_angle(i) / 2)];
                rotation_2 = rotation_1';

                head1_direction = rotation_1 * main_vector_direction(:, i);
                head2_direction = rotation_2 * main_vector_direction(:, i);

                firstHead = [1 i];
                middleHead = [2 i];
                lastHead = [3 i];

                headsX(firstHead(1), firstHead(2)) = x_targ(i) + head1_direction(1) * arrow_side_length(i);
                headsY(firstHead(1), firstHead(2)) = y_targ(i) + head1_direction(2) * arrow_side_length(i);
                headsX(middleHead(1), middleHead(2)) = x_targ(i);
                headsY(middleHead(1), middleHead(2)) = y_targ(i);
                headsX(lastHead(1), lastHead(2)) = x_targ(i) + head2_direction(1) * arrow_side_length(i);
                headsY(lastHead(1), lastHead(2)) = y_targ(i) + head2_direction(2) * arrow_side_length(i);

                firstTail = 3 * (i - 1) + 1;
                lastTail = 3 * (i - 1) + 2;
                tailsX(firstTail) = x_orig(i);
                tailsY(firstTail) = y_orig(i);
                tailsX(lastTail) = x_targ(i);
                tailsY(lastTail) = y_targ(i);
            end
            
            % Plot them
            edges_handle{2} = fill(headsX, headsY, options.edge_color,...% '-k',...
                                   'linewidth', options.edge_thickness,...
                                   'EdgeColor', options.edge_color);
            edges_handle{3} = plot(axis_handle, tailsX, tailsY, '-k',...
                                   'linewidth', options.edge_thickness,...
                                   'Color', options.edge_color);
        end
    end

    %% Nodes
    if numel(options.node_values) == N
%        if grasp_is_octave
%            node_values = fix (size(colormap, 1) / range (node_values) * (node_values - min (node_values)));
%        end
        if size(input_graph.layout, 2) == 2
            nodes_handle = scatter(axis_handle, input_graph.layout(:, 1), input_graph.layout(:, 2), options.node_display_size, options.node_values, 'filled');
            if options.node_marker_edge_width > 0
                set(nodes_handle, 'MarkerEdgeColor', 'k');
                set(nodes_handle, 'LineWidth', options.node_marker_edge_width);
%                nodes_handle.MarkerEdgeColor = 'k';
%                nodes_handle.LineWidth = options.node_marker_edge_width;
            end
        else
            set(gcf,'CurrentAxes', axis_handle)
            nodes_handle = scatter3(input_graph.layout(:, 1), input_graph.layout(:, 2), input_graph.layout(:, 3), options.node_display_size, options.node_values, 'filled');
            if options.node_marker_edge_width > 0
                set(nodes_handle, 'MarkerEdgeColor', 'k');
                set(nodes_handle, 'LineWidth', options.node_marker_edge_width);
%                nodes_handle.MarkerEdgeColor = 'k';
%                nodes_handle.LineWidth = options.node_marker_edge_width;
            end
        end
    else
        if size(input_graph.layout, 2) == 2
            nodes_handle = scatter(axis_handle, input_graph.layout(:, 1), input_graph.layout(:, 2), options.node_display_size, 'b', 'filled');
            if options.node_marker_edge_width > 0
                set(nodes_handle, 'MarkerEdgeColor', 'k');
                set(nodes_handle, 'LineWidth', options.node_marker_edge_width);
%                nodes_handle.MarkerEdgeColor = 'k';
%                nodes_handle.LineWidth = options.node_marker_edge_width;
            end
        else
            set(gcf,'CurrentAxes', axis_handle)
            nodes_handle = scatter3(input_graph.layout(:, 1), input_graph.layout(:, 2), input_graph.layout(:, 3), options.node_display_size, 'b', 'filled');
            if options.node_marker_edge_width > 0
                set(nodes_handle, 'MarkerEdgeColor', 'k');
                set(nodes_handle, 'LineWidth', options.node_marker_edge_width);
%                nodes_handle.MarkerEdgeColor = 'k';
%                nodes_handle.LineWidth = options.node_marker_edge_width;
            end
        end
    end
    
    if numel(options.highlight_nodes) > 0
        scatter(axis_handle,...
                input_graph.layout(options.highlight_nodes, 1), input_graph.layout(options.highlight_nodes, 2),...
                options.highlight_nodes_size, options.highlight_nodes_color, 'o',...
                'LineWidth', options.highlight_nodes_width);
    end
    
    node_text_ID_test = (ischar(options.node_text) && strcmp(options.node_text, 'ID'));
    if numel(options.node_text) > 0 && (numel(options.node_text) ~= grasp_nb_nodes(input_graph) || node_text_ID_test)
        if ~node_text_ID_test && isfield(input_graph, 'node_names') && numel(input_graph.node_names) >= grasp_nb_nodes(input_graph)
            options.node_text = input_graph.node_names;
        else
            if ~node_text_ID_test
                warning('GraSP:showGraphID', 'Using node ID as node_text!');
            end
            options.node_text = arrayfun(@(i) num2str(i), (1:grasp_nb_nodes(input_graph))', 'UniformOutput', false);
        end
    end
    if numel(options.node_text_shift) ~= size(input_graph.layout, 2)
        options.node_text_shift = options.node_text_shift * ones(1, size(input_graph.layout, 2));
    end
    if numel(options.node_text) > 1
        for i = 1:grasp_nb_nodes(input_graph)
            if size(input_graph.layout, 2) == 2
                text(input_graph.layout(i, 1) + options.node_text_shift(1),...
                     input_graph.layout(i, 2) + options.node_text_shift(2),...
                     options.node_text{i},...
                     'FontSize', options.node_text_fontsize,...
                     'BackgroundColor', options.node_text_background_color,...
                     'EdgeColor', options.node_text_background_edge);
            else
                text(input_graph.layout(i, 1) + options.node_text_shift(1),...
                     input_graph.layout(i, 2) + options.node_text_shift(2),...
                     input_graph.layout(i, 3) + options.node_text_shift(3),...
                     options.node_text{i},...
                     'FontSize', options.node_text_fontsize,...
                     'BackgroundColor', options.node_text_background_color,...
                     'EdgeColor', options.node_text_background_edge);
            end
        end
    end
    
    %% Colorbar
    
    if options.show_colorbar
        colorbar(axis_handle);
    else
        colorbar(axis_handle, 'off');
    end
    
    % Show the colorbar for the edge weights => need for a new set of
    % axes
    if options.edge_colorbar && ~isempty(options.edge_colormap)
        axis_edge_cm_handle = axes;
        axis_edge_cm_handle.Title.String = 'edge_colorbar';
        linkprop([axis_handle, axis_edge_cm_handle], {'DataAspectRatio', 'PlotBoxAspectRatio', 'Position', 'XLim', 'YLim', 'ZLim', 'OuterPosition'});

        axis_edge_cm_handle.Visible = 'off';
        axis_edge_cm_handle.XTick = [];
        axis_edge_cm_handle.YTick = [];
        
        colormap(axis_edge_cm_handle, options.edge_colormap);
        caxis(options.edge_color_scale);
        colorbar(axis_edge_cm_handle, 'southoutside');
        axis_edge_cm_handle.Position = axis_handle.Position;
        
        axes(axis_handle);
    end
     
    %% Plot
    set(axis_handle, 'NextPlot', prev_hold);
end