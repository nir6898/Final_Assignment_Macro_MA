%% Macro Theory A - Final Assignment
% Aiyagari Model in Continuous Time with Redistribution
% Main script that generates all required outputs

clear; clc; close all;

% Initialize parameters and set up model
initialize_params;

% Solve the baseline model (τ = 20%)
disp('======= Solving baseline model with τ = 20% =======');
[r_eq, T_eq, Sol] = solve_model(params, grid);

% Display results
disp('======= Model solved =======');
fprintf('Equilibrium interest rate: r = %.4f\n', r_eq);
fprintf('Equilibrium lump-sum transfer: T = %.4f\n', T_eq);
fprintf('Capital demand: K_D = %.4f\n', Sol.K_Dem);
fprintf('Capital supply: K_S = %.4f\n', Sol.Kagg);
fprintf('Aggregate consumption: C = %.4f\n', Sol.Agg_Cons);
fprintf('Aggregate output: Y = %.4f\n', Sol.Agg_output);

% Check capital market clearing and government budget balance
K_error = abs((Sol.K_Dem - Sol.Kagg)/Sol.Kagg);
budget_error = abs((T_eq - params.tau * Sol.Agg_Labor * Sol.wage)/Sol.Agg_output);
fprintf('Capital market clearing error: %.6f\n', K_error);
fprintf('Government budget constraint error: %.6f\n', budget_error);

% Check goods market clearing (Y = C + δK)
goods_market_error = abs((Sol.Agg_output - Sol.Agg_Cons - params.del*Sol.Kagg)/Sol.Agg_output);
fprintf('Goods market clearing error: %.6f\n', goods_market_error);

% Generate required plots
plot_model_results(grid, params, Sol);
