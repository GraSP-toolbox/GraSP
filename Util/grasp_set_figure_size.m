%Set figure size.
%
%When figures are docked, this function resizes all figures.
%
%   GRASP_SET_FIGURE_SIZE(figure_handle, new_size) resizes to new_size =
%       [new_width, new_height] the figure described by figure_handle, or
%       all figures when they are docked. new_size unit is 'pixels'.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2019)
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

function grasp_set_figure_size(figure_handle, new_size)
    tmp_units = get(figure_handle, 'Units');
    set(figure_handle, 'Units', 'pixels');
    cur_pos = get(figure_handle, 'Position');
    if sum(cur_pos(3:4) == new_size) == 2
        return;
    end
    switch get(figure_handle, 'WindowStyle')
        case 'docked'
            desktop   = com.mathworks.mde.desk.MLDesktop.getInstance;
            drawnow;
            container = desktop.getGroupContainer('Figures').getTopLevelAncestor;
            % This is becoming tricky: we set the widget size to the target 
            % figure size, and then adjust the widget size until the figure 
            % size matches the target size (because it will be off).
            widget_dimensions = new_size;
            container.setSize(java.awt.Dimension(widget_dimensions(1), widget_dimensions(2)));
            drawnow
            cur_pos = get(figure_handle, 'Position');
            while cur_pos(3) ~= new_size(1)
                diff = new_size(1) - cur_pos(3);
                delta_width = sign(diff) * ceil(abs(diff / 2));
                widget_dimensions(1) = widget_dimensions(1) + delta_width;
                container.setSize(java.awt.Dimension(widget_dimensions(1), widget_dimensions(2)));
                drawnow
                cur_pos = get(figure_handle, 'Position');
            end
            while cur_pos(4) ~= new_size(2)
                diff = new_size(2) - cur_pos(4);
                delta_height = sign(diff) * ceil(abs(diff / 2));
                widget_dimensions(2) = widget_dimensions(2) + delta_height;
                container.setSize(java.awt.Dimension(widget_dimensions(1), widget_dimensions(2)));
                drawnow
                cur_pos = get(figure_handle, 'Position');
            end
        case 'normal'
            cur_pos = get(figure_handle, 'Position');
            cur_pos(3:4) = new_size;
            set(figure_handle, 'Position', cur_pos);
        otherwise
            error('Unknown ''WindowStyle'' figure property.');
    end
    set(figure_handle, 'Units', tmp_units);
end