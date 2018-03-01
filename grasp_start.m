%Sets Matlab search path to use this toolbox.
%
%   GRASP_START()
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

function grasp_start()
    pwd = [fileparts(mfilename('fullpath')), filesep];
    addpath(pwd);
    
    tmp_dir = '3rdParty';
    addpath([pwd, tmp_dir]);
    dep_list = grasp_dependencies_list;
    rmpath([pwd, tmp_dir]);

    dirs = {'Duality',...
            'Graphs',...
            'Graphs/Tools',...
            'Operators',...
            'Plotting',...
            'Signals',...
            'Stats',...
            'Util'};
    
    for k = 1:numel(dep_list)
        if numel(dep_list(k).optional) > 0 && dep_list(k).optional
            continue;
        end
        dir = dep_list(k).name;
        if numel(dep_list(k).root_dir) > 0
            dir = [dir, dep_list(k).root_dir];
        end
        root_path = [pwd, '3rdParty/', dir];
        for p = 1:numel(dep_list(k).path_list)
            addpath([root_path, dep_list(k).path_list{p}]);
        end
    end
    
    for k = 1:numel(dirs)
        addpath([pwd, dirs{k}]);
    end
    
    global GRASP_OPT_TOOLS
    GRASP_OPT_TOOLS = [];
end