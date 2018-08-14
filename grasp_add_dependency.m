%Add a local dependency. This dependency is necessarily optional and will
%need to be started with GRASP_START_OPT_3RD_PARTY.
%
% Format of the structure is described below (starred fields are
% mandatory):
%          url*: to download the toolbox (zipfile)
%         name*: of the toolbox
%      root_dir: if the toolbox archive has a root folder, use it,
%                otherwise, leave empty (root_dir will be 'name/')
%         debug: set to 1 to disable 
%   init_script: if there is a script to execute to install or initialize
%                the toolbox, place it here (it should be recheable from
%                root_dir or any of the directories listed in path_list)
%      mex_flag: whether init_script requires compilation of mex files
%     path_list: a list of subdirs of root_dir to include in Matlab path
%       ref_bib: a bibliography reference (see GRASP_BIBLIOGRAPHY)
%       patches: cell with a set of couple of strings
%                (file_to_patch_relative to root_dir, patch file in
%                3rdParty) to patch the toolbox after being uncompressed.
%  start_script: if there is a script to execute each time matlab is
%                started, set it here
%       mex_dir: directory where mex files are
%         mexes: cell array of string cells where each string is an input
%                to the function MEX.
%  dependencies: set of dependencies of this dependency (i.e. toolboxes to
%                start alongside this one)
%
%   GRASP_ADD_DEPENDENCY(new_dependency) add new_dependency to the mat file
%       '3rdParty/local_dependencies.mat'.
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

function grasp_add_dependency(new_dependency)
    %% Path to the 3rdParty folder
    pwd = [fileparts(mfilename('fullpath')), filesep];
    thirdparty_dir = [pwd, '3rdParty'];
    
    %% Get the base list
    addpath(thirdparty_dir);
    grasp_dep_list = grasp_dependencies_list;
    rmpath(thirdparty_dir);
    
    %% Remove useless fields
    fn1 = fieldnames(grasp_dep_list(1));
    fn2 = fieldnames(new_dependency);
    tf = ismember(fn2, fn1);
    clean_new_dependency = struct();
    for K = 1:length(tf)
        if tf(K)
            clean_new_dependency.(fn2{K}) = new_dependency.(fn2{K});
        end
    end
    
    %% Force optional
    clean_new_dependency.optional = 1;
    
    %% Cleanup
    if clean_new_dependency.name(end) ~= '/'
        clean_new_dependency.name(end + 1) = '/';
    end
    if sum(clean_new_dependency.name == '/') > 1
        error('The field ''name'' should not contain any slash (''/'')!');
    end
    
    %% Add path
    if ~isfield(clean_new_dependency, 'path_list')
        clean_new_dependency.path_list = {'.'};
    end
    
    %% Save
    local_dep_file = [thirdparty_dir '/local_dependencies.mat'];
    if exist(local_dep_file, 'file')
        load(local_dep_file, 'local_dependencies');
        local_dependencies{end + 1} = clean_new_dependency; %#ok
    else
        local_dependencies{1} = clean_new_dependency;     %#ok 
    end
    save(local_dep_file, 'local_dependencies');
end