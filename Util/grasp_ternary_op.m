%Classical ternary operator.
%
%Note: calling GRASP_TERNARY_OP(test, out_true, out_false) directly
%involves Matlab computing out_true and out_false prior to calling 
%GRASP_TERNARY_OP. This may not be desirable, and a lazier approach needed
%(for example, if either out_true or out_false is undefined, or return
%errors). In which case, it is more efficient to use EVAL:
%EVAL(GRASP_TERNARY_OP(test, 'out_true_code', 'out_false_code'))
%
%   out = GRASP_TERNARY_OP(test, out_true, out_false) ternary operator
%       returninig out_true if test is true, and out_false otherwise.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Southern California, USA
% (2017-2019).
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

function out = grasp_ternary_op(test, out_true, out_false)
    if test
        out = out_true;
    else
        out = out_false;
    end
end