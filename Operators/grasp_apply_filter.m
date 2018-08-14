%Applies a graph filter given a graph.
%
%   H = GRASP_APPLY_FILTER(graph, filter) computes the matrix H of the
%       graph filter filter on graph.
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
                out = polyvalm(filter.data, graph.Z);
            case 'kernel'
                out = grasp_fourier_inverse(graph, diag(filter.data(graph.eigvals)));
            case 'convolution'
                out = grasp_convolution(graph, filter.data);
            case 'matrix'
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