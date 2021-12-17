%Constructs a graph from the Barabási–Albert model. The graph is non
%directed and unweighted.
%
%Note: default GRASP_SHOW_GRAPH settings when options.build_layout is set
%to true show edges of the initial graph in black and new edges in gray.
%However, the graph is still unweighted and default options can be
%overwritten.
%
%   graph = GRASP_BARABASI_ALBERT(N, m) constructs a Barabási–Albert graph
%       with N nodes, initial Erdős-Renyi graph with m < N nodes, wiring
%       probability p = 0.5, and m edges added for each new node.
%
%   GRASP_BARABASI_ALBERT(..., options) optional parameters:
%
%   options.initial_graph: use an initial complete graph with 
%       options.initial_graph <= m nodes. If initial_graph is a grasp graph
%       structure, with m <= nb_nodes, this graph will be used (default: 0
%       for and Erdős-Renyi graph with m nodes).
%
%   options.build_layout: false to keep the layout undefined, 'graphviz'
%       to compute the layout using GRASP_LAYOUT, 'spectral' to compute
%       it using GRASP_LAYOUT_SPECTRAL, true to keep the initial graph
%       layout and arange the new nodes in a circle around it (default:
%       true).
%
%   options.layout_distance_factor: if build_layout is true, defines how
%       far the new nodes are from the initial graph (default: 2).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Southern California, USA
% (2018-2020).
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

function graph = grasp_barabasi_albert(N, m, varargin)
    %% Parameters
    default_param = struct(...
        'initial_graph', 0,...
        'build_layout', true,...
        'layout_distance_factor', 2);
    if isempty(varargin)
        options = struct;
    elseif numel(varargin) > 1
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Initialization
    
    graph = grasp_struct;
    graph.A = zeros(N);
    
    if m < 1 || mod(m, 1) > 0
        error('GraSP:BarabasiAlbert:M', 'm must be at least 1, and integer.');
    end

    % Note: m is the number of edges connecting a new node to the rest of the graph
    % while m0 is the number of nodes in the initial graph.
    if isstruct(options.initial_graph)
        if isfield(options.initial_graph, 'A')
            if grasp_is_directed(options.initial_graph)
                error('GraSP:BarabasiAlbert:InitGraphDirected', 'The intial graph in options.initial_graph should be non directed.');
            end
            g_init = options.initial_graph;
            m0 = grasp_nb_nodes(g_init);
        else
            error('GraSP:BarabasiAlbert:InitGraph', 'options.initial_graph should be an integer > 0 or a valid grasp graph structure.');
        end
    elseif isnumeric(options.initial_graph) && options.initial_graph >= 0 && mod(options.initial_graph, 1) == 0
        if options.initial_graph > 0
            m0 = options.initial_graph;
            g_init = grasp_complete(m0);
        elseif m == 1
            % Only one node is needed (ER with one node is a trivial graph)
            g_init = grasp_struct;
            g_init.A = 0;
            g_init.layout = [5 5];
            m0 = 1;
        else
            % m > 1
            % There is no do/while loop in matlab, hence the trick of the first test below.
            g_init = [];
            % The second test ensures that the graph has no isolated node.
            % FIXME: do we also need to check if the graph is connected?
            while isempty(g_init) || sum(diag(grasp_degrees(g_init)) == zeros(m, 1)) > 0
                g_init = grasp_erdos_renyi(m, 0.5, 'directed', false);
            end
            m0 = m;
        end
    else
        error('GraSP:BarabasiAlbert:InitGraph', 'options.initial_graph should be an integer > 0 or a valid grasp graph structure.');
    end
    
    % m0 must be at least m, otherwise, there is not enough nodes to connect a new node
    % to m of the previous nodes.
    if m0 < m
        error('GraSP:BarabasiAlbert:InitGraphM0', 'The initial graph must have at least m = %d nodes.', m);
    end
    graph.A(1:m0, 1:m0) = g_init.A;
    
    %% Add nodes
    for i = (m0 + 1):N
        %% Probability distribution
        cumul = (blkdiag(tril(ones(i - 1)), zeros(N - i + 1)) * sum(graph.A, 2))';
        cumul = cumul / sum(sum(graph.A));
        
        %% Draw neighbors without replacement
        new_neighbors = zeros(1, m);
        for j = 1:m
            while new_neighbors(j) == 0 || sum(new_neighbors(1:(j-1)) == new_neighbors(j)) > 0
                r = rand(1);
                new_neighbors(j) = 1;
                while cumul(new_neighbors(j)) <= r
                    new_neighbors(j) = new_neighbors(j) + 1;
                end
            end
        end
        
        %% Add neighbors
        for j = new_neighbors
            graph.A(j, i) = 1;
            graph.A(i, j) = 1;
        end
    end
    
    %% Sparsify
    % TODO: build a sparse graph from the beginning for scalability
    graph.A = sparse(graph.A);
    
    %% Layout
    switch options.build_layout
        case false
            return
        case 'graphviz'
            graph.layout = grasp_layout(graph);
        case 'spectral'
            graph.layout = grasp_layout_spectral(graph);
        otherwise
            if size(g_init.layout, 2) > 2
                g_init.layout = g_init.layout(:, 1:2);
            end

            k = N - m0;
            center = mean(g_init.layout);
            r = options.layout_distance_factor * sqrt(max(sum((g_init.layout - center) .^ 2, 2)));
            if r == 0
                r = 4; % Case of m0 = 1, for example.
            end
            z = r * exp(1i * 2 * pi * (0:(k - 1)) / k);
            graph.layout = [g_init.layout ; center + [real(z') imag(z')]];
            
            graph.A_layout = graph.A;
            graph.A_layout(1:m0, 1:m0) = 2 * graph.A_layout(1:m0, 1:m0);
            graph.show_graph_options.layout_boundaries = 0.05;
            graph.show_graph_options.edge_colormap = flipud(colormap('gray'));
            graph.show_graph_options.edge_thickness = 2;
            graph.show_graph_options.edge_color_scale = [0.5 2];
    end
end
