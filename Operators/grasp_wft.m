%Computes the windowed Fourier transform as defined by [Shuman et al.
%2015].
%
%   wft = GRASP_WFT(graph, f) computes the wft matrix such that
%   wft(i, k) = <f,g_{i,k}> with g_{i,k} the atom centered around node i
%   and at frequency lambda_k. The function g used to define the window is
%   a gaussian kernel of decay 3.
%
%   GRASP_WFT(..., kernel_decay) use a gaussian window of prescribed decay.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015-11-01)
% 
% benjamin.girault@ens-lyon.fr
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

function wft = grasp_wft(graph, f, kernel_decay)
    %% Initializations
    N = grasp_nb_nodes(graph);
    wft(N, N) = 0;
    if nargin == 2
        kernel_decay = 3;
    end
    kernel = grasp_heat_kernel(graph, kernel_decay);
    kernel = kernel / norm(kernel);

    %% Windowed Fourier Transform
    for i = 1:N
        window = grasp_generalized_translation(graph, i, kernel);
        for k = 1:N
            atom = grasp_generalized_modulation(graph, k, window);
            wft(i, k) = sum(f .* conj(atom));
        end
    end
end