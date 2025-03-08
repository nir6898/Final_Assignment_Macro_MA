function [r_eq, T_eq, Sol, iter] = solve_model(params, grid)
% SOLVE_MODEL Solves the Aiyagari model with redistribution using Broyden's method
%   This function implements the algorithm described in the assignment:
%   1. Guess r and T
%   2. Solve firm's problem
%   3. Solve household's problem (HJB)
%   4. Find stationary distribution (KFE)
%   5. Check market clearing and government budget
%   6. Update r and T using Broyden's method
%
%   Inputs:
%     params: Struct with model parameters
%     grid:   Struct with grids and matrices
%
%   Outputs:
%     r_eq:   Equilibrium interest rate
%     T_eq:   Equilibrium transfer
%     Sol:    Struct with solution objects (policies, distribution, etc.)

% Set tolerance for convergence
tolerance = 1e-4;
max_iter = 100;

% Initialize Broyden's method
% x = [r; T]
x = [0.03; 0.05]; % Initial guess
x_old = x;
f_old = market_clearing(x, grid, params);

% Identity matrix as initial Jacobian approximation
B = eye(2);

for iter = 1:max_iter
    % Compute market clearing condition
    f = market_clearing(x, grid, params);
    
    % Check for convergence
    if max(abs(f)) < tolerance
        % disp(['Convergence achieved! Iteration: ', num2str(iter)]);
        break;
    end
    
    % Report current state
    % disp(['Iteration: ', num2str(iter), ', r = ', num2str(x(1)), ...
    %       ', T = ', num2str(x(2)), ', Errors = [', num2str(f(1)), ', ', num2str(f(2)), ']']);
    
    % Compute step using Broyden's method
    s = x - x_old;
    y = f - f_old;
    
    % Update Jacobian approximation using Broyden's update
    if norm(s) > 0
        B = B + ((y - B*s) * s') / (s'*s);
    end
    
    % Update variables for next iteration
    x_old = x;
    f_old = f;
    
    % Compute next guess
    x = x - B\f;
    
    % Ensure r is within reasonable bounds
    x(1) = max(min(x(1), params.rho - 0.001), -params.del + 0.001);
    
    % Ensure T is positive
    x(2) = max(x(2), 0.001);
end

% Final equilibrium values
r_eq = x(1);
T_eq = x(2);

% Get full solution at equilibrium
[~, Sol] = model_func(grid, params, r_eq, T_eq);

end

function [f, Sol] = market_clearing(x, grid, params)
% MARKET_CLEARING Computes market clearing errors for given r and T
%   This function solves the model for given r and T and returns:
%   1. Capital market clearing error
%   2. Government budget balance error

r = x(1);
T = x(2);

% Solve the model
[~, Sol] = model_func(grid, params, r, T);

% Calculate relative excess demand for capital
K_error = (Sol.K_Dem - Sol.Kagg) / Sol.Kagg;

% Calculate government budget balance error
gov_income = params.tau * Sol.Agg_Labor * Sol.wage;
budget_error = (T - gov_income) / Sol.Agg_output;

% Return both errors for Broyden's method
f = [K_error; budget_error];

end
