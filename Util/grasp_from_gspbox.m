%Convert a GSPbox graph structure to a GraSP graph structure.
%
%Requires both GraSP and the GSPbox to be started.
%
%   graph = GRASP_FROM_GSPBOX(gsp_graph) convert gsp_graph to a GraSP
%   structure.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2019)
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

function graph = grasp_from_gspbox(gsp_graph)
    if (gsp_graph.directed)
        error('Full support for directed graphs is not finished in GraSP!');
    end
    
    graph = grasp_struct;
    
    graph.A = gsp_graph.W;
    graph.layout = gsp_graph.coords * 5 + 5;
    switch (gsp_graph.lap_type)
        case 'combinatorial'
            graph.fourier_version = 'standard laplacian';
            graph.M = gsp_graph.L;
            graph.Finv = gsp_graph.U;
            graph.F = graph.Finv';
            graph.eigvals = gsp_graph.e;
        case 'normalized'
            graph.fourier_version = 'normalized laplcian';
            graph.M = gsp_graph.L;
            graph.Finv = gsp_graph.U;
            graph.F = graph.Finv';
            graph.eigvals = gsp_graph.e;
    end
    if strcmp(gsp_graph.type, 'GraSP')
        graph.A_layout = gsp_graph.grasp_A_layout;
        graph.distances = gsp_graph.grasp_distances;
        graph.node_names = gsp_graph.grasp_node_names;
        graph.T = gsp_graph.grasp_T;
        graph.background = gsp_graph.grasp_background;
        graph.Q = gsp_graph.grasp_Q;
        graph.Z = gsp_graph.grasp_Z;
    end
    graph.gsp_plotting = gsp_graph.plotting;
end