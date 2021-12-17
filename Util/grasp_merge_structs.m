%Set defaults values to a structure
%
%   s_out = GRASP_MERGE_STRUCTS(s_defaults, s_out) set the default values
%       of s_out to s_defaults.
%
%   ... = GRASP_MERGE_STRUCTS(..., verbose) sets the verbosity (whether or
%       not to show a warning when s_out has a field not in s_defaults)
%       (default: true).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016-2019)
% Copyright Benjamin Girault, École Nationale de la Statistique et de
% l'Analyse de l'Information, Bruz, FRANCE (2020-2021)
% 
% benjamin.girault@ens-lyon.fr
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

function s_out = grasp_merge_structs(s_defaults, s_out, verbose)
    if nargin == 2
        verbose = true;
    end
    fields = fieldnames(s_defaults);
    for i = 1:numel(fields)
        if ~isfield(s_out, fields{i})
            s_out.(fields{i}) = s_defaults.(fields{i});
        end
    end
    if verbose
        cellfun(@(f) warning('GraSP:MergeStruct:MissingDefaultField', ['Field ''' f ''' has no default value, it may be incorrect.']), setxor(fields, fieldnames(s_out)));
    end
end
