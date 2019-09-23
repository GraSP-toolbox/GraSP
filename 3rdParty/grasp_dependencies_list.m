%Necessary information on third party Matlab libraries.
%
%   list = GRASP_DEPENDENCIES_LIST() structure array of third party
%   dependencies with url, name, root_dir, init_script to run for
%   install, and path_list.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016-2019)
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

function list = grasp_dependencies_list
    list = [];
    
    cur_dep = 1;
    
    list(cur_dep).url = 'https://github.com/bgirault-usc/MyPatcher/archive/v1.1.1.zip';
    list(cur_dep).name = 'MyPatcher/';
    list(cur_dep).root_dir = 'MyPatcher-1.1.1/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).debug = 0;
    list(cur_dep).mex_flag = [];
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/4266/versions/4/download/zip';
    list(cur_dep).name = 'GrTheory/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/fileexchange/submissions/10922/v/2/download/zip';
    list(cur_dep).name = 'MatlabBGL/';
    list(cur_dep).root_dir = 'matlab_bgl/';
    list(cur_dep).path_list = {'.', 'graphs'};
    list(cur_dep).optional = 1;
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/10548/versions/1/download/zip';
    list(cur_dep).name = 'anneal/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).patches = {'anneal/anneal.m', 'anneal.patch'};
    list(cur_dep).optional = 1;
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/epfl-lts2/gspbox/releases/download/0.6.0/gspbox-0.6.0.zip';
    list(cur_dep).name = 'gspbox/';
    list(cur_dep).root_dir = 'gspbox/';
    list(cur_dep).init_script = 'gsp_install';
    list(cur_dep).start_script = 'gsp_start';
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = 'https://arxiv.org/abs/1408.5781';
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://wiki.epfl.ch/sgwt/documents/sgwt_toolbox-1.02.zip';
    list(cur_dep).name = 'sgwt/';
    list(cur_dep).root_dir = 'sgwt_toolbox/';
    list(cur_dep).mex_dir = 'mex/';
    list(cur_dep).mexes = {{'-largeArrayDims', 'cheby_op_adjoint_mex.cpp', 'cheby_op.cpp'}, {'-largeArrayDims', 'cheby_op_mex.cpp', 'cheby_op.cpp'}};
    list(cur_dep).start_script = 'sgwt_setpath';
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = 'http://dx.doi.org/10.1016/j.acha.2010.04.005';
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/STAC-USC/Active_SSL_with_Sampling_Theory/archive/39bdfe03d3f54fa1d368958829e6898772c6eba7.zip';
    list(cur_dep).name = 'usc_ssl_sampling/';
    list(cur_dep).root_dir = 'Active_SSL_with_Sampling_Theory-39bdfe03d3f54fa1d368958829e6898772c6eba7/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = 'https://doi.org/10.1145/2623330.2623760';
    list(cur_dep).dependencies = {'sgwt'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/STAC-USC/Graph_Learning/archive/GLL-v2.1.zip';
    list(cur_dep).name = 'usc_graph_learning/';
    list(cur_dep).root_dir = 'Graph_Learning-GLL-v2.1/';
    list(cur_dep).path_list = {'.', 'functions', 'misc'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = {'https://doi.org/10.1109/JSTSP.2017.2726975', 'https://doi.org/10.1109/TSIPN.2018.2872157'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/STAC-USC/Disc-GLasso/archive/7d255cad9f2a73a58af2e92ece07e26358e1d143.zip';
    list(cur_dep).name = 'usc_disc_glasso/';
    list(cur_dep).root_dir = 'Disc-GLasso-7d255cad9f2a73a58af2e92ece07e26358e1d143/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = 'https://doi.org/10.1109/ICASSP.2017.7952698';
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/STAC-USC/symmetric_grid/archive/d3150ee84e370e40ca5868d867f2bf744d846d80.zip';
    list(cur_dep).name = 'usc_symmetric_grid/';
    list(cur_dep).root_dir = 'symmetric_grid-d3150ee84e370e40ca5868d867f2bf744d846d80/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = 'https://doi.org/10.1109/ICASSP.2017.7952929';
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/STAC-USC/GraphStructures/archive/master.zip';
    list(cur_dep).name = 'usc_graphs/';
    list(cur_dep).root_dir = 'GraphStructures-master/';
    list(cur_dep).path_list = {'ToyGraphGSPExample/'};
    list(cur_dep).optional = 1;
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/ychtanaka/FastGSSS/archive/master.zip';
    list(cur_dep).name = 'fast_gsss/';
    list(cur_dep).root_dir = 'FastGSSS-master/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = 'https://doi.org/10.1109/TSP.2019.2908129';
    list(cur_dep).dependencies = {'sgwt'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/STAC-USC/graph_learning_properties/archive/master.zip';
    list(cur_dep).name = 'usc_graph_learning_properties/';
    list(cur_dep).root_dir = 'graph_learning_properties-master/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = {'https://doi.org/10.1109/JSTSP.2017.2726975', 'https://doi.org/10.1109/TSP.2018.2813337'};
    list(cur_dep).dependencies = {'usc_graph_learning'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/bgirault-usc/Molene-Dataset/archive/master.zip';
    list(cur_dep).name = 'bgirault_molene_dataset/';
    list(cur_dep).root_dir = 'Molene-Dataset-master/';
    list(cur_dep).path_list = {'Matlab/'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = {'https://doi.org/10.1109/EUSIPCO.2015.7362637'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'https://github.com/chebfun/chebfun/archive/v5.7.0.zip';
    list(cur_dep).name = 'chebfun/';
    list(cur_dep).root_dir = 'chebfun-5.7.0/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).optional = 1;
    list(cur_dep).ref_bib = {'https://www.chebfun.org/docs/guide/chebfun_guide.pdf'};
end
