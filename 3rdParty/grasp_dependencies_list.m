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
% Angeles, California, USA (2016)
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
    list(cur_dep).debug = 0;
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://perso.ens-lyon.fr/patrick.flandrin/pack_emd.zip';
    list(cur_dep).name = 'EMD/';
    list(cur_dep).root_dir = 'package_emd/';
    list(cur_dep).init_script = 'install_emd';
    list(cur_dep).path_list = {'.', 'EMDs', 'utils'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://github.com/altmany/export_fig/archive/f0af704d84608f5a69e3de82581869e7b6161d4f.zip';
    list(cur_dep).name = 'ExportFig/';
    list(cur_dep).root_dir = 'export_fig-f0af704d84608f5a69e3de82581869e7b6161d4f/';
    list(cur_dep).path_list = {'.'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/fileexchange/submissions/39275/v/4/download/zip';
    list(cur_dep).name = 'HistogramDistance/';
    list(cur_dep).path_list = {'.'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/4266/versions/4/download/zip';
    list(cur_dep).name = 'GrTheory/';
    list(cur_dep).path_list = {'.'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/fileexchange/submissions/10922/v/2/download/zip';
    list(cur_dep).name = 'MathBFL/';
    list(cur_dep).path_list = {'matlab_bgl', 'matlab_bgl/graphs'};
    
    cur_dep = cur_dep + 1;
    
    list(cur_dep).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/10548/versions/1/download/zip';
    list(cur_dep).name = 'anneal/';
    list(cur_dep).path_list = {'.'};
    list(cur_dep).patches = {'anneal/anneal.m', 'anneal.patch'};
    
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
end
