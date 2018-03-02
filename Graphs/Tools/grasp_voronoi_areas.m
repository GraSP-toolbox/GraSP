%Given a graph layout, compute the Voronoi diagram, and the area of each
%cell of the diagram.
%
%   VA = GRASP_VORONOI_AREAS(graph) compute for each vertex the area of
%   its Voronoi cell.
%
%   GRASP_VORONOI_AREAS(..., diagram_ax_handle) show the Voronoi diagram
%   on the axes provided.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Southern California, USA
% (2017-2018).
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

function VA = grasp_voronoi_areas(graph, diagram_ax_handle)
    % Trick: for each vertex, add four vertices mirrored by one of the
    % boundary edges, thus ensuring that each cell is within the boundary
    coords = graph.layout;
    N = size(coords, 1);
    coords(5 * N, 2) = 0;
    for i = 1:N
        coords(N + i, :) = [-coords(i, 1) coords(i, 2)];
        coords(2 * N + i, :) = [coords(i, 1) -coords(i, 2)];
        coords(3 * N + i, :) = [20 - coords(i, 1) coords(i, 2)];
        coords(4 * N + i, :) = [coords(i, 1) 20 - coords(i, 2)];
    end

    % Voronoi Cells
    [voro_verts, voro_cells] = voronoin(coords);
    
    % Plot these cells
    if nargin > 1 && diagram_ax_handle ~= 0
        plot(diagram_ax_handle, graph.layout(:, 1), graph.layout(:, 2), 'xr');
        hold on
        for i = 1:N
            plot(diagram_ax_handle,...
                 [voro_verts(voro_cells{i}, 1) ; voro_verts(voro_cells{i}(1), 1)]',...
                 [voro_verts(voro_cells{i}, 2) ; voro_verts(voro_cells{i}(1), 2)]',...
                 '-b');
        end
        hold off
        xlim(diagram_ax_handle, [-2 12]);
        ylim(diagram_ax_handle, [-2 12]);
    end
    
    % Compute the Voronoi areas
    VA = cellfun(@(c) polyarea(voro_verts(c, 1), voro_verts(c, 2)), voro_cells(1:N));
end