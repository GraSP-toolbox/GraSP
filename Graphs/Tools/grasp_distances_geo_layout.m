%Construct an approximate distance matrix from geographic coordinates. Uses
%the Haversine formula, in combination with an approximate latitude
%dependent earth radius.
%
%   graph = GRASP_DISTANCES_GEO_LAYOUT(graph) computes the approximate
%       distance matrix in kilometers.
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

function distances = grasp_distances_geo_layout(graph)
    R = approximate_radius(graph.layout(:, 2));
    N = grasp_nb_nodes(graph);
    tmp_1 = repmat(1:N, N, 1);
    tmp_2 = tmp_1';
    distances = haversine_approx(graph.layout(tmp_1(:), 2), graph.layout(tmp_1(:), 1), graph.layout(tmp_2(:), 2), graph.layout(tmp_2(:), 1), R);
    distances = reshape(distances, N, N);
end

function [dist_haversine,dist_angular,dist_lat,dist_lon] = haversine_approx(latitude_1, longitude_1, latitude_2, longitude_2, radius)
    dist_lat = deg2rad(latitude_2 - latitude_1);
    dist_lon = deg2rad(longitude_2 - longitude_1);
    latitude_1 = deg2rad(latitude_1);
    latitude_2 = deg2rad(latitude_2);
    dist_angular = sqrt((sin(dist_lat ./ 2)) .^ 2 + cos(latitude_1) .* cos(latitude_2) .* (sin(dist_lon ./ 2)) .^ 2);
    dist_haversine = 2 * radius .* asin(dist_angular);
end

function R = approximate_radius(latitudes)
    R_equator = 6378;
    R_pole = 6357;
    phi = (max(latitudes) + min(latitudes)) / 2;
    s = sin(deg2rad(phi));
    c = cos(deg2rad(phi));
    R = sqrt(((R_equator ^ 2 * c) ^ 2 + (R_pole ^ 2 * s) ^ 2) / ((R_equator * c) ^ 2 + (R_pole * s) ^ 2));
end