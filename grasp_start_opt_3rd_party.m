%Start an optional third party toolbox. Starred entries correspond to local
%dependencies.
%
%   GRASP_START_OPT_3RD_PARTY() lists all optional toolboxes
%
%   GRASP_START_OPT_3RD_PARTY(toolbox_id) starts the given toolbox_id.
%       toolbox_id can be either the id the toolbox or its name (as given
%       by GRASP_START_OPT_3RD_PARTY()).
%
%   opt_tools = GRASP_START_OPT_3RD_PARTY(0) returns opt_tools, the struct
%       array with info on the toolboxes (useful only to 
%       GRASP_BIBLIOGRAPHY).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016-2018)
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

function opt_tools_ret = grasp_start_opt_3rd_party(toolbox_id)
    %% Build the list
    pwd = [fileparts(mfilename('fullpath')), filesep];
    tmp_dir = [pwd, '3rdParty'];
    addpath(tmp_dir);
    dep_list = grasp_dependencies_list;
    rmpath(tmp_dir);
    nb_base_dependencies = numel(dep_list);
    
    % Local dependencies
    mfile_dir = [fileparts(mfilename('fullpath')), filesep];
    local_dep_file = [mfile_dir '/3rdParty/local_dependencies.mat'];
    if exist(local_dep_file, 'file')
        load(local_dep_file, 'local_dependencies');
        for cur_dep = 1:numel(local_dependencies)
            cur_soft_id = numel(dep_list) + 1;
            cur_dep_fields = fieldnames(local_dependencies{cur_dep});
            for f_id = 1:numel(cur_dep_fields)
                dep_list(cur_soft_id).(cur_dep_fields{f_id}) = local_dependencies{cur_dep}.(cur_dep_fields{f_id});
            end
        end
    end
    
    % Cleanup
    opt_tools = repmat(dep_list(1), sum([dep_list.optional]), 1);
    cur_tool = 1;
    for i = 1:numel(dep_list)
        if numel(dep_list(i).optional) > 0 && dep_list(i).optional
            opt_tools(cur_tool).name = dep_list(i).name(1:(end - 1)); % removing the trailing /
            opt_tools(cur_tool).id = cur_tool;
            opt_tools(cur_tool).dep_id = i;
            opt_tools(cur_tool).url = dep_list(i).ref_bib;
            opt_tools(cur_tool).root_dir = [];
            cur_tool = cur_tool + 1;
        end
    end

    %% Do we list or start?
    if nargin == 0
        for i = 1:numel(opt_tools)
            if opt_tools(i).dep_id <= nb_base_dependencies
                fprintf('#%d\t%s\n', opt_tools(i).id, opt_tools(i).name);
            else
                fprintf('#%d*\t%s\n', opt_tools(i).id, opt_tools(i).name);
            end
        end
        return
    end
    
    %% toolbox_id = 0: return the array
    if toolbox_id == 0
        opt_tools_ret = opt_tools;
        return;
    end
    
    %% Starting then...
    global GRASP_OPT_TOOLS
    if numel(GRASP_OPT_TOOLS) ~= numel(opt_tools)
        GRASP_OPT_TOOLS = zeros(numel(opt_tools), 1);
    end
    if ischar(toolbox_id)
        actual_id = -1;
        for i = 1:numel(opt_tools)
            if strcmpi(opt_tools(i).name, toolbox_id)
                actual_id = i;
                break;
            end
        end
        toolbox_id = actual_id;
    end
    if toolbox_id < 0 || toolbox_id > numel(opt_tools)
        error('No such optional toolbox!');
    end
    if GRASP_OPT_TOOLS(toolbox_id)
        % Already started
        return;
    end

    soft = dep_list(opt_tools(toolbox_id).dep_id);
    
    for d = 1:numel(soft.dependencies)
        grasp_start_opt_3rd_party(soft.dependencies{d})
    end
    
    root_path = [pwd, '3rdParty/', soft.name, soft.root_dir];
    if numel(soft.start_script) > 0
        addpath(root_path);
        start = str2func(soft.start_script);
        start();
    end
    for p = 1:numel(soft.path_list)
        addpath([root_path, soft.path_list{p}]);
    end
    GRASP_OPT_TOOLS(toolbox_id) = 1;
    
    disp(['Toolbox ''' opt_tools(toolbox_id).name ''' started.']);
end
