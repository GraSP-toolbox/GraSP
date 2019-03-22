%Generate a gif image file from a graph and several graph signals.
%
%   GRASP_GENERATE_GIF(fh, filename, graph, signals, titles, title_font_size, show_graph_options)
%       generates the GIF file filename by iteratively plotting 
%       signals(:, i) using graph. titles{i} is used to give a title to the 
%       figure, with font size title_font_size. The figure handle fh is 
%       used.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018)
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

function grasp_generate_gif(fh, filename, graph, signals, titles, title_font_size, show_graph_options, varargin)
    %% Parameters
    default_param = struct(...
        'signals_delays', 0.1,...
        'overlay_function', [],...
        'overlay_update', []);
    if nargin == 7
        options = struct;
    elseif nargin > 8
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Initializations
    ah = get(fh, 'CurrentAxes');
%     disp('<show_graph>');
    nh = grasp_show_graph(ah, graph, show_graph_options);
%     disp('</show_graph>');
    overlay_h = [];
    if numel(options.overlay_function) > 0
        hold(ah, 'on');
        overlay_h = options.overlay_function(ah);
        hold(ah, 'off');
    end
    colorbar;
%         disp('<set>');
    set(fh, 'color', 'w');
%         disp('</set>');
    if numel(options.signals_delays) == 1
        options.signals_delays = options.signals_delays * ones(size(signals, 2), 1);
    end
    
    %% Code
%         disp('<set>');
    set(nh, 'CData', signals(:, 1));
%         disp('</set>');
    if numel(options.overlay_update) > 0
        options.overlay_update(overlay_h, 1);
    end
    drawnow;
    f = getframe(fh);
    if ~grasp_is_octave()
        [im, map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        [im, map] = rgb2ind(f.cdata);
    end
    imwrite(im, map, filename, 'DelayTime', options.signals_delays(1), 'LoopCount', inf);
    for k = 2:size(signals, 2)
        title(ah, titles{k}, 'interpreter', 'latex', 'fontsize', title_font_size);
%         disp('<set>');
        set(nh, 'CData', signals(:, k));
%         disp('</set>');
        if numel(options.overlay_update) > 0
            options.overlay_update(overlay_h, k);
        end
        drawnow;
        f = getframe(fh);
        if ~grasp_is_octave()
            [im, map] = rgb2ind(f.cdata, 256, 'nodither');
        else
            [im, map] = rgb2ind(f.cdata);
        end
        imwrite(im, map, filename, 'DelayTime', options.signals_delays(k), 'WriteMode', 'append');
    end
end
