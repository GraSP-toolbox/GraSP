%Computes the graph translation [Girault et al. 2015, IEEE SPL].
%
%   T = GRASP_TRANSLATION(graph) uses the graph structure to build its
%   graph translation.
%
%   T = GRASP_TRANSLATION(..., freqs) use the given frequencies instead of
%   the ones from [Girault et al. 2015, IEEE SPL]. The resulting operator
%   verifies T\chi_l=exp(-i * 2 * pi * freqs(l)) \chi_l.
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

function T = grasp_translation(graph, freqs)
    %% Use the provided freqs if given
    if nargin == 2
        T = grasp_fourier_inverse(graph, diag(exp(-1i * 2 * pi * freqs)));
        return;
    end
    
    %% Intialization
    if ~strcmp(graph.fourier_version, 'standard laplacian') && ~strcmp(graph.fourier_version, 'normalized laplacian')
        error('Error: The Fourier transform should be the result of the decomposition of a Laplacian matrix!');
    end
    
    %% Eigenvalue upper bound
    rho = 2;
    if strcmp(graph.fourier_version, 'standard laplacian')
        rho = grasp_lapl_eigval_upper_bound(graph);
    end
    
    %% Translation operator
    T = expm(-1i * pi * sqrtm(full(graph.L) / rho));
end