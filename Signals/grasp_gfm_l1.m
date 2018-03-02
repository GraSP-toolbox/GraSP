%Compute the graph Fourier modes with an l1-norm graph variation.
%
%   U = GRASP_GFM_L1(graph, M, Q, K) compute the K first (Delta,Q)-graph
%   Fourier modes, with Delta(x)=||Mx||_1.
%
%   [..., fvals, eflags, outputs, lambdas] GRASP_GFM_L1(...) returns also
%   the result of the optimization functions (FMINCON).
%
% Authors:
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, University of Southern California, USA
% (2017-2018).
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

function [U, fvals, eflags, outputs, lambdas] = grasp_gfm_l1(graph, M, Q, K)
    %% Initializations
    N       = grasp_nb_nodes(graph);
    
    U       = zeros(N, K);
    fvals   = cell(K, 1);
    eflags  = cell(K, 1);
    outputs = cell(K, 1);
    lambdas = cell(K, 1);

    %% Compute the GFMs one by one
    for cur_mode = 1:K
        %% Build the linear program with quadratic constraints
        % solve for y = [u' x']'
        %   x our target graph signal
        %   u >= |Mx|
        % and minimize 1' * u
    
        Qobj = sparse(2 * N, 2 * N);
        fobj = [ones(N, 1) ; zeros(N, 1)];
        cobj = 0;

        H   = cell(2 * N, 1);
        k   = cell(2 * N, 1);
        d   = cell(2 * N, 1);
        Heq = cell(cur_mode - 1 + 1, 1);
        keq = cell(cur_mode - 1 + 1, 1);
        deq = cell(cur_mode - 1 + 1, 1);

        % Constraint on u for u >= |Mx|
        K1 = [-eye(N)  M];  % K1 * y <= 0  <=> -u + Mx <= 0  <=>  u >=   Mx
        K2 = [-eye(N) -M];  % K2 * y <= 0  <=> -u - Mx <= 0  <=>  u >= - Mx
        for jj = 1:N
            H{jj} = sparse(2 * N, 2 * N);
            k{jj} = K1(jj, :)';
            d{jj} = 0;

            H{jj + N} = sparse(2 * N, 2 * N);
            k{jj + N} = K2(jj, :)';
            d{jj + N} = 0;
        end

        % Constraint on x for x orthogonal to the previous modes
        for prev_mode = 1:(cur_mode - 1)
            Heq{prev_mode} = sparse(2 * N, 2 * N);
            keq{prev_mode} = [zeros(N, 1) ; Q * U(:, prev_mode)];
            deq{prev_mode} = 0;
        end

        % Constraint for a unit norm vector
        Heq{cur_mode} = 2 * sparse(blkdiag(zeros(N), Q));
        keq{cur_mode} = zeros(2 * N, 1);
        deq{cur_mode} = -1;


        % Initial vector (random x, projected on the orthogonal space, 
        % unit norm, with u = |x| + epsilon)
        x0 = rand(N, 1);
        for jj = 1:(cur_mode - 1)
            x0 = x0 - U(:, jj)' * Q * x0 * U(:, jj);
        end
        x0 = x0 / sqrt(x0' * Q * x0);
        u0 = abs(M * x0) + 1e-4;
        y0 = [u0 ; x0];


        %% Optim

        options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
                                         'SpecifyObjectiveGradient', true,...
                                         'SpecifyConstraintGradient', true,...
                                         'ConstraintTolerance', 1e-15,...
                                         'OptimalityTolerance', 1e-10,...
                                         'FiniteDifferenceType', 'central',...
                                         'MaxFunctionEvaluations', 5000,...
                                         'Display', 'off',...
                                         'HessianFcn', @(x,lambda) quadhess(x, lambda, Qobj, H, Heq));


        fun = @(x) quadobj(x, Qobj, fobj, cobj);
        nonlconstr = @(x) quadconstr(x, H, k, d, Heq, keq, deq);

        [cstr0, cstreq0] = nonlconstr(y0);
        if sum(cstr0 > 0) + sum(abs(cstreq0) > 1e-2) > 0
            error('infeasible initial point!');
        end


        [y, fvals{cur_mode}, eflags{cur_mode}, outputs{cur_mode}, lambdas{cur_mode}]...
                = fmincon(fun, y0, [], [], [], [], [], [], nonlconstr, options);
        
        if eflags{cur_mode} ~= 1
            warning('fmincon did not converge to a feasible point!');
        end
        
        %% Retrive the solution
        U(:, cur_mode) = y(N + (1:N));
    end

end

%% Helper functions (for a generic quadratic problem with quadratic constraints)

% Objective function, with gradient
function [y, grady] = quadobj(x, Q, f, c)
    y = 1/2 * x' * Q * x + f' * x + c;
    if nargout > 1
        grady = Q * x + f;
    end
end

% Inequality and equality constraints, with gradients
function [y, yeq, grady, gradyeq] = quadconstr(x, H, k, d, Heq, keq, deq)
    nb_ineq_constr = length(H);
    
    nb_eq_constr = 0;
    if nargin > 5
        nb_eq_constr = length(Heq); % jj is the number of inequality constraints
    end
    
    y = zeros(1, nb_ineq_constr);
    for i = 1:nb_ineq_constr
        y(i) = 1/2 * x' * H{i} * x + k{i}' * x + d{i};
    end
    
    yeq = zeros(1, nb_eq_constr);
    for i = 1:nb_eq_constr
        yeq(i) = 1/2 * x' * Heq{i} * x + keq{i}' * x + deq{i};
    end

    if nargout > 2
        grady = zeros(length(x), nb_ineq_constr);
        for i = 1:nb_ineq_constr
            grady(:, i) = H{i} * x + k{i};
        end
        
        gradyeq = zeros(length(x), nb_eq_constr);
        for i = 1:nb_eq_constr
            gradyeq(:, i) = Heq{i} * x + keq{i};
        end
    end
end

% Hessian
function hess = quadhess(~, lambda, Q, H, Heq)
    hess = Q;
    nb_ineq_constr = length(H);
    for i = 1:nb_ineq_constr
        hess = hess + lambda.ineqnonlin(i) * H{i};
    end
    nb_eq_constr = length(Heq);
    for i = 1:nb_eq_constr
        hess = hess + lambda.eqnonlin(i) * Heq{i};
    end
end