%Convert a GraSP graph structure to a GSPbox graph structure.
%
%Requires both GraSP and the GSPbox to be started.
%
%   gsp_graph = GRASP_TO_GSPBOX(graph) convert graph to a GSPbox structure.
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

function gsp_graph = grasp_to_gspbox(graph)
    gsp_graph.W = graph.A;
    gsp_graph.grasp_A_layout = graph.A_layout;
    gsp_graph.coords = (graph.layout - 5) / 5;
    gsp_graph.grasp_distances = graph.distances;
    gsp_graph.grasp_node_names = graph.node_names;
    gsp_graph.L = graph.L;
    switch (graph.fourier_version)
        case 'standard laplacian'
            gsp_graph.lap_type = 'combinatorial';
        case 'normalized laplacian'
            gsp_graph.lap_type = 'normalized';
        case 'adja'
            error('Unsupported Fourier tranform by GSPbox!');
    end
    if numel(graph.eigvals) > 1
        gsp_graph.e = graph.eigvals;
        gsp_graph.U = graph.Finv;
    end
    gsp_graph.grasp_T = graph.T;
    gsp_graph.grasp_background = graph.background;
    gsp_graph.type = 'GraSP';
    
    gsp_graph = gsp_graph_default_parameters(gsp_graph);
end