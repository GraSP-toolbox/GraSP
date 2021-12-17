%Transforms a graph filter of type 'kernel' into a polynomial graph filter
%using polynomial interpolation.
%
%   filter = GRASP_FILTER_KERNEL_TO_POLY(filter) returns filter
%       interpolated into a polynomial filter of degree 5, with samples
%       equidistributed in [0, 2]
%
%   GRASP_FILTER_KERNEL_TO_POLY(..., options) optional parameters:
%
%   options.algorithm: Use one of the following algorithm to obtain a 
%       polynomial:
%       - 'regular_interpolation' (default): polynomial interpolation with
%       equidistributed samples in options.freqs_bounds.
%       - 'graph_fit_lsqr': use the graph frequencies instead of
%       equidistributed samples and perform a least squares approximation
%       - 'graph_interpolation': Lagrange polynomial based on the graph
%       frequencies
%       - 'taylor': use a Taylor expansion around options.lambda0.
%       - 'puiseux': use a Puiseux expansion around options.lambda0.
%   
%   options.poly_degree: defines the degree of the polynomial (default: 5)
%
%   options.graph: for options.type = 'graph_interpolation', the graph
%       structure defining the graph frequencies (default: empty).
%
%   options.freqs_bounds: the interval to compute the polynomial
%       interpolation on (default: [0 2]). See options.type = 
%       'regular_interpolation'.
%
%   options.lambda0: to perform an expansion (default: 0). See options.type
%       = 'taylor' or 'puiseux'.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2019)
% Copyright Benjamin Girault, École Nationale de la Statistique et de
% l'Analyse de l'Information, Bruz, FRANCE (2020-2021)
% 
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

function filter = grasp_filter_kernel_to_poly(filter, varargin)
    %% Parameters
    default_param = struct(...
        'poly_degree', -1,...
        'algorithm', 'regular_interpolation',...
        'graph', [],...
        'freqs_bounds', [0 2],...
        'lambda0', 0);
    if isempty(varargin)
        options = struct;
    elseif length(varargin) > 1
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Checks
    if ~strcmp(filter.type, 'kernel')
        error('filter needs to be a kernelized filter!');
    end
    
    if options.poly_degree == -1
        if ~strcmp(options.algorithm, 'graph_interpoloation')
            options.poly_degree = 5;
        end
    end
    
    %% Core
    filter.type = 'polynomial';
    
    switch options.algorithm
        case 'regular_interpolation'
            freqs = (0:options.poly_degree) / options.poly_degree * (options.freqs_bounds(2) - options.freqs_bounds(1)) + options.freqs_bounds(1);
            filter.data = polyfit(freqs, filter.data(freqs), options.poly_degree);
        case 'graph_interpolation'
            if isempty(options.graph)
                error('Missing ''graph'' parameter!');
            end
            if ~isfield(options.graph, 'eigvals') || numel(options.graph.eigvals) ~= grasp_nb_nodes(options.graph)
                error('Missing graph frequencies! Use grasp_eigendecomposition first!');
            end
            if options.poly_degree == -1
                options.poly_degree =  grasp_nb_nodes(options.graph);
            end
            N = grasp_nb_nodes(options.graph);
            x = options.graph.eigvals;
            y = filter.data(options.graph.eigvals);
            filter.data = zeros(1, N);
            for l = 1:N
                cur_x = x([1:(l - 1) (l + 1):N]);
                filter.data = filter.data + y(l) * prod(x(l) - cur_x) * poly(cur_x);
            end
        case 'graph_fit_lsqr'
            if isempty(options.graph)
                error('Missing ''graph'' parameter!');
            end
            if ~isfield(options.graph, 'eigvals') || numel(options.graph.eigvals) ~= grasp_nb_nodes(options.graph)
                error('Missing graph frequencies! Use grasp_eigendecomposition first!');
            end
            if options.poly_degree == -1
                options.poly_degree =  grasp_nb_nodes(options.graph);
            end
            filter.data = polyfit(options.graph.eigvals, filter.data(options.graph.eigvals), options.poly_degree);
        case 'taylor'
            syms x;
            filter.data = sym2poly(taylor(filter.data(x), x, 'Order', options.poly_degree + 1, 'ExpansionPoint', options.lambda0));
        case 'puiseux'
            syms x;
            filter.data = sym2poly(series(filter.data(x), x, 'Order', options.poly_degree + 1, 'ExpansionPoint', options.lambda0));
        otherwise
            error(['Unknown interpolation algorithm ''' options.algorithm '''!']);
    end
end
