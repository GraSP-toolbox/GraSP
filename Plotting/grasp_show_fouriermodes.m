%Shows some eigenvectors on a graph
%
%   GRASP_SHOW_FOURIERMODES(graph) show on separate subfigures the Fourier
%   modes on the graph.
%
%   GRASP_SHOW_FOURIERMODES(..., options) optional parameters:
%
%   options.modes: displays only the modes indices of the array modes.
%
%   options.titles: titles to the subgraphs. titles is a cell array, that
%       must be of the size the number of modes, or 0 to display the
%       eigenvalues. Any other size will disable titles.
%
%   options.titles_fontsize: sets the font size for the GFT mode titles
%       (default: Matlab default).
%
%   options.titles_margin: Percentage of the whole figure vertical size
%       above an axis that is reserved for the title (default: 0.02).
%
%   options.eigval_precision: when using options.titles = 0, controls how
%       many floating point digits are displayed (default: 2).
%
%   options.energy: shows the energy of the modes.
%
%   options.value_scale: sets the color scale of the plot.
%
%   options.show_edges: whether or not edges are shown on the plots.
%
%   options.cmap: colormap to use.
%
%   options.node_size: node size (scatter).
%
%   options.*: any of the options to GRASP_SHOW_GRAPH.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
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

function grasp_show_fouriermodes(graph, varargin)
    %% Parameters
    default_param = struct(...
        'modes', 0,...
        'titles', cell(1),...
        'titles_fontsize', get(0, 'DefaultAxesTitleFontSize'),...
        'titles_margin', 0.02,...
        'eigval_precision', 2,...
        'energy', 0,...
        'value_scale', [-1 1],...
        'show_edges', true,...
        'cmap', 'default',...
        'node_size', 1000);
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Initializations
    N = grasp_nb_nodes(graph);
    
    if numel(options.modes) == 1 && options.modes == 0
        options.modes = 1:N;
    end
    
    nb_modes = numel(options.modes);
    
    [sa_nbl, sa_nbc] = grasp_subaxis_matrix_dimensions(gcf, nb_modes, graph.background);
    
    if ~isfield(graph, 'layout') || numel(graph.layout) <= 1
        graph.layout = grasp_layout(graph);
    end
    
    if numel(options.titles) == 1
        options.titles = cell(nb_modes, 1);
        for i = 1:nb_modes
            options.titles{i} = sprintf(sprintf('\\\\lambda_{%%d}=%%.%df', options.eigval_precision), options.modes(i), graph.eigvals(options.modes(i)));
        end
    end
    
    if numel(options.value_scale) ~= 2
        m = max(max(abs(graph.Finv(:, options.modes))));
        options.value_scale = [-m m];
    end
    
    %% Graphs
    fh = gcf;
    clf(fh);
    
    options_gsg = rmfield(options, fieldnames(default_param));
    options_gsg.node_display_size = options.node_size;
    options_gsg.value_scale = options.value_scale;
    options_gsg.show_edges = options.show_edges;
    
    for mode_id = 1:nb_modes
        mode = options.modes(mode_id);
        if grasp_is_octave
            curAxis = subplot(sa_nbl, sa_nbc, mode_id);
        else
            if numel(options.titles) ~= nb_modes
                curAxis = grasp_subaxis(fh, sa_nbl, sa_nbc, mode_id, 'S', 0.02, 'P', 0, 'M', 0.02);
            else
                curAxis = grasp_subaxis(fh, sa_nbl, sa_nbc, mode_id, 'SV', 0.01, 'SH', 0, 'P', 0.01, 'M', 0, 'MT', 0.02, 'PT', options.titles_margin);
            end
        end
        
        options_gsg.node_values = graph.Finv(:, mode);
        
        if options.energy
            grasp_show_graph(...
                curAxis,...
                graph,...
                grasp_merge_structs(options_gsg,...
                struct(...
                    'node_values', abs(graph.Finv(:, mode)) .^ 2,...
                    'value_scale', [0 max(abs(options.value_scale))])));
%                     'node_display_size', options.node_size,...
%                     'node_values', abs(graph.Finv(:, mode)) .^ 2,...
%                     'value_scale', [0 max(abs(options.value_scale))],...
%                     'show_edges', options.show_edges)));
        else
            grasp_show_graph(...
                curAxis,...
                graph,...
                options_gsg);
%                 grasp_merge_structs(options_gsg,...
%                 struct(...
%                     'node_display_size', options.node_size,...
%                     'node_values', graph.Finv(:, mode),...
%                     'value_scale', options.value_scale,...
%                     'show_edges', options.show_edges)));
        end
        if strcmp(options.cmap, 'default')
            if ~options.energy && options.value_scale(1) == -options.value_scale(2)
                grasp_set_blue_red_colormap(curAxis);
            end
        else
            colormap(curAxis, options.cmap);
        end
        if nargin >= 2 && numel(options.titles) == nb_modes
            title(curAxis, options.titles{mode_id}, 'FontSize', options.titles_fontsize);
        end
    end
end