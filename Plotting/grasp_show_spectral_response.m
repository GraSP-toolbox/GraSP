%Plots a graph filter spectral response, and its interpolation if
%available.
%
%   GRASP_SHOW_SPECTRAL_RESPONSE(axis_handle, graph, filter) plots the
%       frequency response and, if available its interpolation in the
%       spectral domain.
%
%   GRASP_SHOW_SPECTRAL_RESPONSE(..., signal) also plot a graph signal
%       spectral components before and after filtering.
%
%   [freq_resp_handle, legend_handle, freq_resp_interp_handle] = GRASP_SHOW_SPECTRAL_RESPONSE(...)
%       Also return the graphical elements handle to tweak them with no
%       input signal (these are respectively the frequency response, the 
%       input signal, the output signal, the legend and the interpolated 
%       frequency response).
%
%   [freq_resp_handle, legend_handle, input_handle, output_handle, freq_resp_interp_handle] = GRASP_SHOW_SPECTRAL_RESPONSE(...)
%       Also return the graphical elements handle to tweak them when a
%       signal is provided as input to the function
%       (these are respectively the frequency response, the input signal,
%       the output signal, the legend and the interpolated frequency
%       response).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018)
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

function [freq_resp_handle, legend_handle, varargout] = grasp_show_spectral_response(axis_handle, graph, filter, signal)
    %% Checks
    if numel(graph.eigvals) ~= grasp_nb_nodes(graph)
        error('Missing graph Fourier transform! please use grasp_eigendecomposition first.');
    end
    
    %% Start holding plots
    prev_hold = get(axis_handle, 'NextPlot');
    if strcmp(prev_hold, 'replace')
        cla(axis_handle, 'reset');
    end
    hold(axis_handle, 'on');
    
    %% Spectral response
    x = 0:(graph.eigvals(end) / 10000):graph.eigvals(end);
    switch filter.type
        case 'polynomial'
            freq_resp_interp_handle = plot(axis_handle, x, polyval(filter.data, x), '-');
            freq_resp_handle        = scatter(axis_handle, graph.eigvals, polyval(filter.data, graph.eigvals), 'b');
            
            varargout{2 * double(nargin == 4) + 1} = freq_resp_interp_handle;
        case 'chebpoly'
            freq_resp_interp_handle = plot(axis_handle, x, grasp_apply_filter(x, filter), '-');
            freq_resp_handle        = scatter(axis_handle, graph.eigvals, grasp_apply_filter(graph.eigvals, filter), 'b');
            
            varargout{2 * double(nargin == 4) + 1} = freq_resp_interp_handle;
        case 'kernel'
            freq_resp_interp_handle = plot(axis_handle, x, filter.data(x), '-');
            freq_resp_handle        = scatter(axis_handle, graph.eigvals, filter.data(graph.eigvals), 'b');
            
            varargout{2 * double(nargin == 4) + 1} = freq_resp_interp_handle;
        case 'convolution'
            freq_resp_handle        = plot(axis_handle, graph.eigvals, grasp_fourier(graph, filter.data), ':o');
        case 'matrix'
            error('Only graph filters with a spectral response are supported!');
        otherwise
            error(['Unknown filter type ''' filter.type '''!']);
    end
    
    %% Input / Output signals
    if nargin == 4
        input  = grasp_fourier(graph, signal);
        output = grasp_fourier(graph, grasp_apply_filter(graph, filter, signal));

        yyaxis(axis_handle, 'right')
        input_handle  = stem(graph.eigvals, input, ':');
        output_handle = stem(graph.eigvals, output, 'x');
        legend_handle = legend([freq_resp_handle, input_handle, output_handle], {'Frequency response', 'Input', 'Output'});

        tmp = ylim;
        tmp = abs(tmp(1) / (tmp(2) - tmp(1)));

        yyaxis(axis_handle, 'left');
        ylim(axis_handle, [-tmp / (1 - tmp) 1]);
        
        varargout{1} = input_handle;
        varargout{2} = output_handle;
    else
        legend_handle = legend(freq_resp_handle, {'Frequency response'});
    end
    
    %% Stop holding plots if necessary
    set(axis_handle, 'NextPlot', prev_hold);
end
