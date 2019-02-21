%Builds a graph signal with prescribed spectral correlation matrix and mean
%and normal distribution on each node.
%
%   signals = GRASP_INTRINSIC_STATIONARY_REALIZATION(graph, variogram, sill) 
%       builds a single realization of an intrinsic stationary random field
%       with prescribed variogram that is a function of distance (and its
%       limit when the distance goes to infinity, the sill). The field is 
%       zero-mean.
%
%   signals = GRASP_INTRINSIC_STATIONARY_REALIZATION(..., mean) same with a
%       prescribed mean.
%
%   signals = GRASP_INTRINSIC_STATIONARY_REALIZATION(..., nb_realizations) 
%       same with multiple realizations.
%
%   [..., Cov] = GRASP_INTRINSIC_STATIONARY_REALIZATION(...) also returns
%       the true covariance matrix of the random graph signal.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2018-2019)
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

function [signals, Cov] = grasp_intrinsic_stationary_realization(graph, variogram, sill, mean, nb_realizations)
    %% Preprocessing
    N = grasp_nb_nodes(graph);
    
    if nargin == 3
        mean = zeros(N, 1);
        nb_realizations = 1;
    elseif nargin == 4
        if numel(mean) ~= N
            nb_realizations = mean(1);
            mean = zeros(N, 1);
        else
            nb_realizations = 1;
        end
    end
    if numel(mean) ~= N
        error('GraSP:IntrinsicStationaryRealization:MeanWrongSize', '''options.mean'' should be of size N, the number of vertices.');
    end
    
    %% Covariance
    Cov = sill - arrayfun(variogram, graph.distances);
    
    %% White noise
    white_noise = normrnd(0, 1, N, nb_realizations);
    
    %% White noise transform matrix to get the right covariance matrix
    try
        L = chol(Cov, 'lower');
    catch
        warning('off', 'MATLAB:sqrtm:SingularMatrix');
        A = sqrtm(Cov);
        if sum(sum(imag(A))) > 0
            error('GraSP:IntrinsicStationaryRealization:NegativeMatrix', 'Negative matrix.');
        end
        [~, R] = qr(A);
        L = R';
    end
    
    %%
    signals = mean + L * white_noise;
end