%Applies a graph filter given a graph.
%
%   H = GRASP_APPLY_FILTER(graph, filter) computes the matrix H of the
%       graph filter filter on graph.
%
%   y = GRASP_APPLY_FILTER(x, filter) applies filter to each entry of x
%       independently (not for 'convolution' and 'matrix' filter types).
%
%   output = GRASP_APPLY_FILTER(..., input) computes the output of the
%       filter on the graph signal input.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018)
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

function out = grasp_apply_filter(graph, filter, signal)
    if nargin == 2
        switch filter.type
            case 'polynomial'
                if isstruct(graph)
                    out = polyvalm(filter.data, graph.Z);
                else
                    out = polyval(filter.data, graph);
                end
            case 'chebpoly'
                map_shift = (filter.data.interval(2) + filter.data.interval(1)) / 2;
                map_scale = 2 / (filter.data.interval(2) - filter.data.interval(1));
                
                if isstruct(graph)
                    U0 = speye(grasp_nb_nodes(graph));
                    tmp = map_scale * (graph.Z - map_shift * U0);
                else
                    U0 = ones(size(graph));
                    tmp = map_scale * (graph - map_shift * U0);
                end
                U1 = (filter.data.cheb_kind) * tmp;
                out = filter.data.coeffs(1) * U0;
                
                if numel(filter.data.coeffs) > 1
                    out = out + filter.data.coeffs(2) * U1;
                    for m = 3:numel(filter.data.coeffs)
                        if isstruct(graph)
                            new_U = 2 * tmp * U1 - U0;
                        else
                            new_U = 2 * tmp .* U1 - U0;
                        end
                        U0 = U1;
                        U1 = new_U;
                        out = out + filter.data.coeffs(m) * U1;
                    end
                end
            case 'kernel'
                if isstruct(graph)
                    out = grasp_fourier_inverse(graph, diag(filter.data(graph.eigvals)));
                else
                    out = filter.data(graph);
                end
            case 'convolution'
                if ~isstruct(graph)
                    error('Graph only filter!');
                end
                out = grasp_convolution(graph, filter.data);
            case 'matrix'
                if ~isstruct(graph)
                    error('Graph only filter!');
                end
                out = filter.data;
            otherwise
                error(['Unknown filter type ''' filter.type '''!']);
        end
    elseif nargin == 3
        switch filter.type
            case 'polynomial'
                tmp = signal;
                out = zeros(numel(signal), 1);
                for i = numel(filter.data):-1:1
                    out = out + filter.data(i) * tmp;
                    if i ~= 1
                        tmp = graph.Z * tmp;
                    end
                end
            case 'chebpoly'
                map_shift = (filter.data.interval(2) + filter.data.interval(1)) / 2;
                map_scale = 2 / (filter.data.interval(2) - filter.data.interval(1));
                
                U0 = signal;
                U1 = (filter.data.cheb_kind) * map_scale * (graph.Z * signal - map_shift * signal);
                out = filter.data.coeffs(1) * U0;
                
                if numel(filter.data.coeffs) > 1
                    out = out + filter.data.coeffs(2) * U1;
                    for m = 3:numel(filter.data.coeffs)
                        new_U = 2 * map_scale * (graph.Z * U1 - map_shift * U1) - U0;
                        U0 = U1;
                        U1 = new_U;
                        out = out + filter.data.coeffs(m) * U1;
                    end
                end
            case 'kernel'
                out = grasp_fourier_inverse(graph, filter.data(graph.eigvals) .* grasp_fourier(graph, signal));
            case 'convolution'
                out = grasp_convolution(graph, filter.data, signal);
            case 'matrix'
                out = filter.data * signal;
            otherwise
                error(['Unknown filter type ''' filter.type '''!']);
        end
    else
        error('Wrong number of arguments!');
    end
end