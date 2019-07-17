%Setup third party Matlab libraries.
%
%   GRASP_INIT_3RD_PARTY(with_local_deps) downloads and install in this 
%       directory third party Matlab libraries.
%
%   GRASP_INIT_3RD_PARTY(options) optional parameters:
%
%   options.with_local_deps: true if libraries added using 
%       GRASP_ADD_DEPENDENCY are also installed. (default: false)
%
%   options.without_mex: disable building of dependencies that require 
%       compilation of mex files (default: false)
%
%   options.reinstall: whether already installed libraries should be
%       reinstalled (default: false)
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016-2018)
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

function grasp_init_3rd_party(varargin)
    %% Parameters
    default_param = struct(...
        'with_local_deps', false,...
        'without_mex', false,...
        'reinstall', false);
    if nargin == 0
        options = struct;
    elseif nargin > 1
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);

    %% List
    softwares = grasp_dependencies_list;

    mfile_dir = [fileparts(mfilename('fullpath')), filesep];
    local_dep_file = [mfile_dir '/local_dependencies.mat'];
    if ~exist(local_dep_file, 'file')
        local_dependencies = [];
        save(local_dep_file, 'local_dependencies');
    end
    if options.with_local_deps
        load(local_dep_file, 'local_dependencies');
        for cur_dep = 1:numel(local_dependencies)                 %#ok
            cur_soft_id = numel(softwares) + 1;
            cur_dep_fields = fields(local_dependencies{cur_dep}); %#ok
            for f_id = 1:numel(cur_dep_fields)
                softwares(cur_soft_id).(cur_dep_fields{f_id}) = local_dependencies{cur_dep}.(cur_dep_fields{f_id});
            end
        end
    end
    
    %% Matlab version
    package_save = @(url, file) urlwrite(url, file);  %#ok
    matlab_ver = ver('Matlab');
    if numel(matlab_ver) > 0 && str2double(matlab_ver.Version) >= 8.4
        package_save = @(url, file) websave(file, url);
    end
    
    %% MyPatcher path
    path_mypatcher = [fileparts(mfilename('fullpath')), filesep, softwares(1).name, softwares(1).root_dir];

    %% Let's start...
    pwd = [fileparts(mfilename('fullpath')), filesep];
    for soft_id = 1:numel(softwares)
        soft = softwares(soft_id);
        if numel(soft.debug) > 0 && soft.debug
            continue;
        end
        if options.without_mex && ((numel(soft.mex_flag) > 0 && soft.mex_flag) || (isfield(soft, 'mexes') && numel(soft.mexes) > 0))
            continue;
        end
        dir = [pwd, soft.name];
        if exist(dir, 'dir')
            if ~options.reinstall
                disp([soft.name(1:(end - 1)) ' already installed, skipping...']);
                continue;
            else
                % delete(dir); % We don't actually remove anything, just in case
            end
        end
        disp(['Installing ' soft.name(1:(end - 1)) '...']);
        package_save(soft.url, 'tmp.zip');
        mkdir(pwd, soft.name);
        unzip('tmp.zip', dir);
        delete('tmp.zip');
        if numel(soft.patches) ~= 0
            addpath(path_mypatcher);
            for i = 1:2:numel(soft.patches)
                mypatcher(soft.patches{i}, soft.patches{i + 1}, soft.patches{i});
            end
            rmpath(path_mypatcher);
        end
        if numel(soft.init_script) ~= 0
            old_path = path;
            toolpath = [fileparts(mfilename('fullpath')), filesep, soft.name, soft.root_dir];
            addpath(toolpath);
            init = str2func(soft.init_script);
            init();
            path(old_path);
        end
        if isfield(soft, 'mexes') && numel(soft.mexes) ~= 0
            curdir = pwd;
            cd([dir soft.root_dir soft.mex_dir]);
            for i = 1:numel(soft.mexes)
                eval(['mex(' sprintf('''%s'', ', soft.mexes{i}{1:end - 1}) '''' soft.mexes{i}{end} ''');'])
            end
            cd(curdir);
        end
    end
end
