%Constructs a complete graph
%
%   graph = GRASP_COMPLETE(N) constructs a complete graph with N nodes.
%
%   options.build_layout: true if the graph should include a layout
%       (default: true). Nodes will be aranged on a circle.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, University of Southern California, USA
% (2020).
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

function graph = grasp_complete(N, varargin)
    %% Parameters
    default_param = struct(...
        'build_layout', true);
    if isempty(varargin)
        options = struct;
    elseif length(varargin) > 1
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Build the graph
    graph = grasp_struct;
    graph.A = ones(N) - eye(N);
    
    if options.build_layout
        z = 0.9 * exp(1i * 2 * pi * (0:(N - 1)) / N);
        graph.layout = [real(z') imag(z')];
        graph.layout = 10 * (graph.layout + 1) / 2;
    end
end
