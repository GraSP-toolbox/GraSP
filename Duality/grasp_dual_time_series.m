%Maps a graph to a time series using various approaches.
%
%   ts = GRASP_DUAL_TIME_SERIES(graph) creates a random time_series ts
%   using [Girault et al. 2014, IEEE ICASSP]
%
%   [ts, P, bins_center, bins_bounds, final_loss, ssl_centroids] = GRASP_DUAL_TIME_SERIES(...)
%   returns also the random walk matrix P, the center of bins for
%   amplitudes bins_center and their bounds bins_bounds. final_loss is the
%   cost value of the simulated annealing solution, and ssl_centroids are
%   the two centroids used for GSSL.
%
%   GRASP_DUAL_TIME_SERIES(..., options)  optional parameters:
%
%   options.nb_timesteps: number of samples of the time series (default:
%       1000).
%
%   options.method: either 'ssl' [Girault et al. 2014, IEEE ICASSP] or
%       'anneal' [Campanharo et al. 2011, PLOS], or 'custom' (with custom
%       bins definitions). (default: 'ssl').
%
%       GSSL method
%
%   options.adjacency_matrix: type of adjacency matrix (see
%       GRASP_FARTHEST). (default: 'distance').
%
%   options.sigma_ssl: value of paramter sigma
%
%   options.alpha_ssl: value of paramter alpha
%
%       Simulated annealing method
%
%   options.init_temp: inital temperature (default: 20).
%
%   options.stop_temp: stopping temperature (default: 1).
%
%       Custom method
%
%   options.bins_center: bins centers for amplitudes
%
%   options.bins_bounds: bins bounds for amplitudes
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

function [ts, P, bins_center, bins_bounds, final_loss, ssl_centroids] = grasp_dual_time_series(graph, varargin)
    %% Parameters
    default_param = struct(...
        'method', 'ssl',...
        'adjacency_matrix', 'distance',...
     ...% anneal options
        'init_temp', 20,...
        'stop_temp', 1,...
     ...% SSL options
        'sigma_ssl', 1,...
        'alpha_ssl', 2,...
     ...% custom options
        'bins_center', 0,...
        'bins_bounds', 0,...
     ...% time series options
        'nb_timesteps', 1000);
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    
    %% Bins: output
    bins_center = options.bins_center;
    bins_bounds = options.bins_bounds;
    
    %% Check
    if full(sum(graph.A(:) < 0)) > 0
        error('The edge weights must be positive!');
    end
    
    %% Random walk matrix
    D = grasp_degrees(graph);
    P = sparse(diag(diag(D) .^ -1) * graph.A);
    
    %% Amplitude mapping
    N = grasp_nb_nodes(graph);
    switch options.method
        case 'ssl'
            % Symmetrize the graph
            if grasp_is_directed(graph)
                graph.A = (graph.A + graph.A') / 2;
            end
            
            % Centroids
            [i, j] = grasp_farthest(graph, options);
            ssl_centroids = [i j];
            H = grasp_semi_supervised(graph, options.sigma_ssl, 2 / options.alpha_ssl);
                % (alpha = 2 / mu)
            
            % Input
            X(N, 1) = 0;
            X(i) = 1;
            X(j) = -1;
            
            % SSL
            Y_num = H * X;
            Y_denom = H * abs(X);
            
            % bins centers
            bins_center = Y_num ./ Y_denom;
            final_loss = -1;
        case 'anneal'
            % Load the toolbox
            grasp_start_opt_3rd_party('anneal');
            
            % First, we convert the adjacency matrix to a random walk
            % matrix
            graph.A = P;
            
            % Next: Simulated Annealing
            loss = @(p) full(sum(sum(graph.A .* abs(ones(N, 1) * p - p' * ones(1, N)))));
            firstSol = 1:N;
            options_anneal = anneal();
            options_anneal.Generator = @(p, T) generator(loss, N^2, p, T);
            options_anneal.Verbosity = 2;
            options_anneal.MaxSuccess = 20;
            options_anneal.InitTemp = options.init_temp;
            options_anneal.StopTemp = options.stop_temp;
            labels = anneal(loss, firstSol, options_anneal);
            final_loss = loss(labels);
            
            % bins centers
            bins_center = labels' / (N + 1);
        case 'custom'
            if numel(bins_center) ~= N
                error('Please provide correct bin centers for the amplitude mapping!');
            end
        otherwise
            error('Unrecognize method!');
    end
    if numel(bins_bounds) ~= 2 * N
        [~, IX] = sort(bins_center);
        bins_bounds = grasp_duality_bins_bounds(...
            bins_center,...
            bins_center(IX(1)) - (bins_center(IX(2)) - bins_center(IX(1))) / 2,...
            bins_center(IX(N)) + (bins_center(IX(N)) - bins_center(IX(N - 1))) / 2);
    end
    
    %% Cumulative probability matrix
    P_cumul = full(P);
    for j = 2:N
        P_cumul(:, j) = P_cumul(:, j - 1) + P_cumul(:, j);
    end
    
    %% Time series
    ts(options.nb_timesteps) = 0;
    cur_node = randi(grasp_nb_nodes(graph), 1);
    for t = 1:options.nb_timesteps
        ts(t) = bins_bounds(cur_node, 1) + (bins_bounds(cur_node, 2) - bins_bounds(cur_node, 1)) * rand(1);
        next_rand = rand(1);
        cur_node = find(P_cumul(cur_node, :) > next_rand, 1, 'first');
    end
end

function outPerm = generator(loss, nbIter, prevPerm, T)
    %% Intializations
    N = numel(prevPerm);
    newPerm = prevPerm;
    
    %% Gaussian distribution for segment size
    pd = ProbDistUnivParam('Normal', [0 sqrt(T + N)]);
    
    %% Metropolis sampling
    outPerm = prevPerm;
    for i = 1:nbIter
        %% New permutation
        curPerm = outPerm;

        first = randi(N, 1);

        segmentSize = max(1, ceil(abs(icdf(pd, rand(1)))));
        if first + segmentSize - 1 > N
            continue;
        end

        newPosition = randi(N - segmentSize + 1, 1);
        
        last = first + segmentSize - 1;
        newLast = newPosition + segmentSize - 1;
        
        segment = curPerm(first:last);
        tmpPerm = curPerm(setdiff(1:N, first:last));
        newPerm(1:newPosition - 1) = tmpPerm(1:newPosition - 1);
        newPerm(newPosition:newLast) = segment;
        newPerm(newLast + 1:N) = tmpPerm(newPosition:end);
        
        %% Acceptation
        p = exp((loss(curPerm) - loss(newPerm)) / T);
        if rand(1) <= p
            outPerm = newPerm;
        end
    end
end