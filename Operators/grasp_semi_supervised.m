%Computes a graph SSL kernel. [Avrachenkov et al., 2012: Generalized
%Optimization Framework for Graph-based Semi-supervised Learning]
%
%Compared to [Girault et al., 2014:  Semi-Supervised Learning for Graph to
%Signal Mapping: a Graph Signal Wiener Filter Interpretation], alpha=2\mu.
%
%NOTE: This function is faster if the normalized Laplacian
%eigendecomposition has been performed.
%
%   H_SSL = GRASP_SEMI_SUPERVISED(g, sigma, mu) computes the SSL kernel of
%   the graph g with parameters sigma and mu.
%
%   [H_SSL, H_SSL_hat] = GRASP_SEMI_SUPERVISED(...) also returns the
%   operator in the spectral domain.
%
%   output = GRASP_SEMI_SUPERVISED(..., input) applies SSL to input
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

function [H_SSL, H_SSL_hat] = grasp_semi_supervised(g, sigma, mu, input)    
    %% Rescaling (sigma)
    N = grasp_nb_nodes(g);
    D = grasp_degrees(g);
    Dsigma =    sparse(diag(  diag(D) .^ (0.5 - sigma)      ));
    Dsigmainv = sparse(diag(  diag(D) .^ (sigma - 0.5)      ));
    
    %% Input / Output
    if nargin == 4
        input = Dsigmainv * input;
        output = (eye(N) + 2 / mu * grasp_laplacian_normalized(g)) \ input;
        output = Dsigma * output;
        H_SSL = output;
        return;
    end

    %% Initialization Fourier
    if ~isfield(g, 'fourier_version') || ~strcmp(g.fourier_version, 'normalized laplacian')
        g = grasp_eigendecomposition(g, 'matrix', 'norm_lapl');
    end
    
    %% Classifier
    H_SSL_hat = sparse(diag(  mu ./ (mu + 2 * g.eigvals)    ));

    H_SSL = Dsigma * g.Finv * H_SSL_hat * g.F * Dsigmainv;
end