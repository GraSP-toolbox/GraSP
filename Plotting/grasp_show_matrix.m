%Plots a matrix, with either regular square elements, or irregular ones.
%IMSHOW is a faster alternative but requires extra work to scale the plot.
%
%   GRASP_SHOW_MATRIX(axis_handle, M) Plots the matrix elements in the
%   plane.
%
%   GRASP_SHOW_MATRIX(..., X, Y) uses X as the horizontal coordinates of
%   the vertical lines, and Y as the vertical coordinates of teh horizontal
%   lines (there are one more coordinate for each than there are rows /
%   columns).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% 
% benjamin.girault@ens-lyon.fr
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

function grasp_show_matrix(axis_handle, M, X, Y)
    %% Initializations
    nbcols = size(M, 2);
    nblines = size(M, 1);
    if nargin <= 2
        X = 0:nbcols;
        Y = 0:nblines;
    end
    
    polys_x(4, nbcols * nblines) = 0;
    polys_y(4, nbcols * nblines) = 0;
    
    %% Indices of the elements
    [ix, iy] = meshgrid(1:nbcols, 1:nblines);
    xs = ix(:);
    ys = iy(:);
    indices = (xs' - 1) * nblines + ys';
    
    %% Polygones
    polys_x(1, :) = X(xs);
    polys_y(1, :) = Y(ys);
    polys_x(2, :) = X(xs + 1);
    polys_y(2, :) = Y(ys);
    polys_x(3, :) = X(xs + 1);
    polys_y(3, :) = Y(ys + 1);
    polys_x(4, :) = X(xs);
    polys_y(4, :) = Y(ys + 1);
    
    %% Plotting
    axes(axis_handle);
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);
    set(axis_handle, 'ydir', 'reverse');
    set(axis_handle, 'DataAspectRatio', [nblines nbcols 1]);
    hold(axis_handle, 'on');
    patch(polys_x, polys_y, M(indices), 'EdgeColor', 'none');
    hold(axis_handle, 'off');
end