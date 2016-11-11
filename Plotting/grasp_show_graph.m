%Plots a graph and a graph signal. See the options below for details on all
%the capabilities.
%
% Arrow definition
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
%   pointed by axis_handle.
%
%   GRASP_SHOW_GRAPH(..., options) optional parameters:
%
%   options.background: plot over a background image.
%
%   options.node_values: plot the values on the nodes using colors
%       (default: no values).
%
%   options.value_scale: use the color scale provided (default: [min max]
%       of node_values).
%
%   options.highlight_nodes: draw a bigger circle around the given nodes.
%
%   options.color_map: use the provide color map (default: 'default').
%
%   options.node_display_size: use the size provided for the nodes
%       (default: 1500).
%
%   options.node_text: use the provided cell to plot also a label
%       associated to each node (default: cell(0)).
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
%   options.edge_color_scale: together with edge_colormap, set the
%       weight boundaries in the colormap.
%
%   options.edge_thickness: use the thickness provided to draw edges or
%       directed edges (default: 0.5).
%
%   options.tail_back_quantity: use the provided quantity to substract from
%       a directed edge tip to draw an arrow (default: 0.2 for a layout
%       with coordinates in [0 10] without maximizing the window ; use 0.08
%       when maximizing).
%
%   options.head_proportion: use the provided peedge_color_scalercentage to
%       compute the arrow's head length (default: 10%).
%
%   [nodes_handle, edges_handle] = GRASP_SHOW_GRAPH(...) returns the handle
%   to the nodes graphics objects (scatter), and edge graphics object. For
%   the edges, it returns an array of 3 elements: the non directed edges
%   (plot), the arrow heads (fill), and the arrow tails (plot)
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016)
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

function [nodes_handle, edges_handle] = grasp_show_graph(axis_handle, graph, varargin)
    %% Parameters
    default_param = struct(...
        'background', graph.background,...
        'node_values', 0,...
        'value_scale', 0,...
        'highlight_nodes', [],...
        'color_map', 'default',...
        'node_display_size', 1500,...
        'node_text', [],...
        'show_edges', true,...
        'edge_color', [0 0 0],...
        'edge_colormap', '',...
        'edge_color_scale', 0,...
        'edge_thickness', 0.5,...
        'tail_back_quantity', 0.2,...
        'head_proportion', 0.1);
    if nargin == 2
        options = struct;
    elseif nargin > 3
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Intialization
    N = grasp_nb_nodes(graph);
    A_layout = graph.A;
    if isfield(graph, 'A_layout') && numel(graph.A_layout) == N ^ 2
        A_layout = graph.A_layout;
    end
    cla(axis_handle);
    hold(axis_handle, 'on');
    % Setting boundaries now for optimisation
    axis(axis_handle, [0 10 0 10]);
    set(axis_handle, 'XTick', [], 'YTick', []);
    if size(graph.layout, 2) == 3  % 3D layout
        zlim(axis_handle, [0 10]);
        view(axis_handle, [20 45]);
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
    colormap(axis_handle, options.color_map);
    
    %% Background
    if ~isempty(options.background)
        [img, map] = imread(options.background);
        info = imfinfo(options.background);
        if ~isempty(map)
            img = ind2rgb(img, map);
        end
        img(:, :, 1) = mean(img, 3);
        img(:, :, 2) = img(:, :, 1);
        img(:, :, 3) = img(:, :, 1);
        for i = 1:size(img, 3)
            img(:, :, i) = flipud(img(:, :, i));
        end
        imagesc([0 10], [0 10], img);
        set(axis_handle, 'DataAspectRatio', [info.Height info.Width 1]);
    end

    %% Layout
    if ~isfield(graph, 'layout') || size(graph.layout, 1) ~= N;
        graph.layout = grasp_layout(graph);
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
        if size(graph.layout, 2) == 2
            edges_handle{1} = plot(axis_handle,...
                                  [graph.layout(rows', 1)' ; graph.layout(cols', 1)'],...
                                  [graph.layout(rows', 2)' ; graph.layout(cols', 2)'],...
                                  '-',...
                                  'linewidth', options.edge_thickness,...
                                  'Color', options.edge_color);
        else
            set(gcf,'CurrentAxes',axis_handle)
            edges_handle{1} = plot3([graph.layout(rows', 1)' ; graph.layout(cols', 1)'],...
                                    [graph.layout(rows', 2)' ; graph.layout(cols', 2)'],...
                                    [graph.layout(rows', 3)' ; graph.layout(cols', 3)'],...
                                    '-',...
                                    'linewidth', options.edge_thickness,...
                                    'Color', options.edge_color);
        end
        % Adjust edge color if necessary
        if ~isempty(options.edge_colormap)
            w = A_layout(sub2ind(size(A_layout), rows, cols));
            min_w = options.edge_color_scale(1);
            span_w = options.edge_color_scale(2) - options.edge_color_scale(1);
            colors = (w - min_w) / (span_w);
            colors = 1 + (1 - min(max(colors, 0), 1)) * (size(options.edge_colormap, 1) - 1);
            min(floor(colors))
            colors = options.edge_colormap(floor(colors), :);
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
            % Arrows
            x_orig = graph.layout(rows', 1)';
            y_orig = graph.layout(rows', 2)';
            x_targ = graph.layout(cols', 1)';
            y_targ = graph.layout(cols', 2)';
            x = x_targ - x_orig;
            y = y_targ - y_orig;
            length_tail = sqrt(x .^ 2 + y .^ 2) * 10;
            x_targ = x_targ - x ./ length_tail * options.tail_back_quantity;
            y_targ = y_targ - y ./ length_tail * options.tail_back_quantity;
            shortest_length = min(length_tail);

            M = numel(x_orig);
            headsX = NaN * zeros(3, M);
            headsY = NaN * zeros(3, M);
        %     headsX = NaN * zeros(1, 4 * M);
        %     headsY = NaN * zeros(1, 4 * M);
            tailsX = NaN * zeros(1, 3 * M);
            tailsY = NaN * zeros(1, 3 * M);
            for i = 1:M
                alpha = atan(shortest_length / length_tail(i));
                rotation_1 = [cos(alpha / 2) -sin(alpha / 2) ; sin(alpha / 2) cos(alpha / 2)];
                rotation_2 = rotation_1';

                main_vector = -[x_targ(i) - x_orig(i) ; y_targ(i) - y_orig(i)];
                head1 = rotation_1 * main_vector;
                head2 = rotation_2 * main_vector;

                firstHead = [1 i];
                middleHead = [2 i];
                lastHead = [3 i];
        %         firstHead = [1 4 * (i - 1) + 1];
        %         middleHead = [1 4 * (i - 1) + 2];
        %         lastHead = [1 4 * (i - 1) + 3];
        %         % separator = [1 4 * (i - 1) + 4];

                headsX(firstHead(1), firstHead(2)) = x_targ(i) + head1(1) * options.head_proportion;
                headsY(firstHead(1), firstHead(2)) = y_targ(i) + head1(2) * options.head_proportion;
                headsX(middleHead(1), middleHead(2)) = x_targ(i);
                headsY(middleHead(1), middleHead(2)) = y_targ(i);
                headsX(lastHead(1), lastHead(2)) = x_targ(i) + head2(1) * options.head_proportion;
                headsY(lastHead(1), lastHead(2)) = y_targ(i) + head2(2) * options.head_proportion;

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
    if grasp_is_octave
        options.node_display_size = options.node_display_size / 100;
    end

    if numel(options.node_values) == N
%        if grasp_is_octave
%            node_values = fix (size(colormap, 1) / range (node_values) * (node_values - min (node_values)));
%        end
        if size(graph.layout, 2) == 2
            nodes_handle = scatter(axis_handle, graph.layout(:, 1), graph.layout(:, 2), options.node_display_size, options.node_values, '.');
        else
            set(gcf,'CurrentAxes', axis_handle)
            nodes_handle = scatter3(graph.layout(:, 1), graph.layout(:, 2), graph.layout(:, 3), options.node_display_size, options.node_values, '.');
        end
    else
        if size(graph.layout, 2) == 2
            nodes_handle = scatter(axis_handle, graph.layout(:, 1), graph.layout(:, 2), options.node_display_size, 'b', '.');
        else
            set(gcf,'CurrentAxes', axis_handle)
            nodes_handle = scatter3(graph.layout(:, 1), graph.layout(:, 2), graph.layout(:, 3), options.node_display_size, 'b', '.');
        end
    end
    
    if numel(options.highlight_nodes) > 0
        scatter(axis_handle, graph.layout(options.highlight_nodes, 1), graph.layout(options.highlight_nodes, 2), options.node_display_size, 'r', 'o');
    end
    
    if numel(options.node_text) > 0 && numel(options.node_text) ~= grasp_nb_nodes(graph)
        if isfield(graph, 'node_names')
            options.node_text = graph.node_names;
        else
            warning('Field ''node_names'' not set in the graph!');
            options.node_text = default_param.node_text;
        end
    end
    if numel(options.node_text) > 1
        for i = 1:grasp_nb_nodes(graph)
            text(graph.layout(i, 1), graph.layout(i, 2), options.node_text{i});
        end
    end
     
     %% Plot
    hold(axis_handle, 'off');
end