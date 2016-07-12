%Builds a graph signal with prescribed spectral correlation matrix and mean
%and normal distribution on each node.
%
%   signal = GRASP_GAUSSIAN_SIGNAL_REALIZATION(graph, mean, S) using the 
%   spectral correlation matrix S, constructs a realization of the graph
%   signal normal variables. The signal has mean as its mean vector.
%
%   signal = GRASP_GAUSSIAN_SIGNAL_REALIZATION(..., m) generate m
%   realizations.
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

function signal = grasp_gaussian_signal_realization(graph, mean, S, m)
    %% Preprocessing
    N = grasp_nb_nodes(graph);
    mean_hat = grasp_fourier(graph, mean);
    Cov_hat = S - mean_hat * mean_hat';  % Covariance
    
    if nargin == 3
        m = 1;
    end
    
    %% X
    X_hat = normrnd(0, 1, N, m);
    
    %% Y = LX
    try
        L = chol(Cov_hat, 'lower');
    catch
        warning('off', 'MATLAB:sqrtm:SingularMatrix');
        A = sqrtm(Cov_hat);
        if sum(sum(imag(A))) > 0
            error('Negative matrix!');
        end
        [~, R] = qr(A);
        L = R';
    end
    Y_hat = L * X_hat;
    
    %% Z = Y + mean
    signal = graph.Finv * Y_hat + kron(ones(1, m), mean);
end