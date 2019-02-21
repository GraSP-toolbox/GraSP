%Computes the reduced frequencies associated to each Fourier mode of a
%graph. See [Girault et al. 2015, IEEE SPL].
%
%   freqs = GRASP_FREQUENCIES(graph) returns the graph frequencies.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
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

function freqs = grasp_frequencies(graph)
    if ~isfield(graph, 'fourier_version') || strcmp(graph.fourier_version, 'n.a.')
        error('Fourier transform not performed!');
    end
    
    if strcmp(graph.fourier_version, 'standard laplacian') || strcmp(graph.fourier_version, 'irregularity-aware')
        freqs = 0.5 * abs(sqrt(graph.eigvals / grasp_lapl_eigval_upper_bound(graph)));
    elseif strcmp(graph.fourier_version, 'normalized laplacian')
        freqs = 0.5 * abs(sqrt(graph.eigvals / 2));
%     elseif 
%         freqs = 0.5 * abs(sqrt(graph.eigvals / max(graph.eigvals)));
    elseif strcmp(graph.fourier_version, 'graph shift') || strcmp(graph.fourier_version, 'diagonalization')
        freqs = 0.25 * (1 - graph.eigvals / max(abs(graph.eigvals)));
    else
        error('Unknown Fourier transform!');
    end
end