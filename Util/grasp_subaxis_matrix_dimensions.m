%Spreads elements in an array such that the number of empty cells is
%minimized.
%
%   [L, C] = GRASP_SUBAXIS_MATRIX_DIMENSIONS(nb_cells) returns the number
%   of lines and the number of columns to display nb_cells elements. gcf is
%   used to get the size of the current figure window.
%
%   [L, C] = GRASP_SUBAXIS_MATRIX_DIMENSIONS(figure_handle, nb_cells) idem
%   using the provided figure handle.
%
%   GRASP_SUBAXIS_MATRIX_DIMENSIONS(..., background) uses a background
%   image to obtain the relation between width and height instead of the
%   window.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Copyright Benjamin Girault, University of Southern California, USA
% (2016-2018).
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


function [L, C] = grasp_subaxis_matrix_dimensions(figure_handle, nb_cells, background)
    if nargin == 1
        nb_cells = figure_handle;
        figure_handle = gcf;
        background = '';
    elseif nargin == 2
        if ischar(nb_cells)
            background = nb_cells;
            nb_cells = figure_handle;
            figure_handle = gcf;
        else
            background = '';
        end
    end

    if strcmp(background, '')
        pos = get(figure_handle, 'Position');
        width = pos(3);
        height = pos(4);
    else
        info = imfinfo(background);
        width = info.Width;
        height = info.Height;
    end
    edge_length = sqrt(width * height / nb_cells);
    L = floor(height / edge_length);
    C = floor(width / edge_length);
    if L * C >= nb_cells
        return;
    end
    if C >= L
        if C * L + L >= nb_cells
            C = C + 1;
            return;
        end
    end
    if C * L + C >= nb_cells
        L = L + 1;
        return;
    end
    C = C + 1;
    L = L + 1;
    return;
end