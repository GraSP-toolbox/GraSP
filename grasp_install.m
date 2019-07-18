%Download and setup necessary third party tools into the directory
%3rdParty/.
%
%   GRASP_INSTALL() install the base third party toolboxes
%
%   GRASP_INSTALL(options) optional parameters:
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
% Inria, FRANCE (2015-2016)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018-2019)
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

function grasp_install(varargin)
    %% Checking we are in the right folder
    mfile_dir = [fileparts(mfilename('fullpath')), filesep];
    if ~exist([mfile_dir '/3rdParty'], 'dir')
        error('The current folder needs to be the root folder grasp of the GraSP toolbox!');
    end

    %% Adding path to Util folder (for grasp_merge_struct)
    addpath([mfile_dir '/Util']);

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
    
    %% Start install
    prev = cd([mfile_dir '/3rdParty']);
    grasp_init_3rd_party(options);
    cd(prev);
end