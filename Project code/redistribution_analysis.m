%% Redistribution Analysis
% This script implements Section 3 of the assignment by solving
% the model for different tax rates and analyzing the results

clear; clc; close all;

% Initialize parameters and set up model
initialize_params;

% Define tax rates to analyze
tax_rates = [0.05:0.05:0.9];  % 5%, 10%, ..., 90%
num_rates = length(tax_rates);

% Initialize arrays to store results
results.r = zeros(num_rates, 1);          % Interest rate
results.T = zeros(num_rates, 1);          % Transfer
results.Y = zeros(num_rates, 1);          % Output
results.K = zeros(num_rates, 1);          % Capital
results.C = zeros(num_rates, 1);          % Consumption
results.TY_ratio = zeros(num_rates, 1);   % Transfer to output ratio
results.tax_revenue = zeros(num_rates, 1); % Tax revenue
results.cons_var = zeros(num_rates, 1);   % Consumption variance

% Store the original tax rate
original_tau = params.tau;

% Loop through tax rates
disp('======= Starting redistribution analysis =======');
for i = 1:num_rates
    % Set new tax rate
    params.tau = tax_rates(i);
    disp(['Solving model for τ = ', num2str(params.tau*100), '%']);
    
    % Solve the model
    [r_eq, T_eq, Sol, iter_needed] = solve_model(params, grid);
    
    % Store results
    results.r(i) = Sol.net_return;
    results.T(i) = Sol.T;
    results.Y(i) = Sol.Agg_output;
    results.K(i) = Sol.Kagg;
    results.C(i) = Sol.Agg_Cons;
    results.TY_ratio(i) = Sol.T / Sol.Agg_output;
    results.tax_revenue(i) = params.tau * Sol.Agg_Labor * Sol.wage;
    
    % Calculate consumption variance
    % First, compute mean consumption
    mean_c = Sol.Agg_Cons;
    
    % Then compute variance
    var_c = 0;
    for z = 1:params.znum
        inds = z_index_range(z, grid.knum);
        var_c = var_c + sum(Sol.h(inds) .* (Sol.c(inds) - mean_c).^2);
    end
    
    results.cons_var(i) = var_c;
    
    disp(['Completed τ = ', num2str(params.tau*100), '%, r = ', num2str(r_eq, '%.4f'), ...
          ', T = ', num2str(T_eq, '%.4f'), ', Y = ', num2str(Sol.Agg_output, '%.4f')]);
    disp(['Converged in ',num2str(iter_needed), ' itetarions']);

end
% Restore original tax rate
params.tau = original_tau;

%% Plot results
figure('Position', [100, 100, 1200, 800]);

% Create a 3x3 subplot
subplot(3, 3, 1);
plot(tax_rates*100, results.Y, 'o-', 'LineWidth', 2);
title('(a) Aggregate Output (Y)');
xlabel('Tax Rate (%)');
ylabel('Output');
grid on;

subplot(3, 3, 2);
plot(tax_rates*100, results.K, 'o-', 'LineWidth', 2);
title('(b) Aggregate Capital (K)');
xlabel('Tax Rate (%)');
ylabel('Capital');
grid on;

subplot(3, 3, 3);
plot(tax_rates*100, results.C, 'o-', 'LineWidth', 2);
title('(c) Aggregate Consumption (C)');
xlabel('Tax Rate (%)');
ylabel('Consumption');
grid on;

subplot(3, 3, 4);
plot(tax_rates*100, results.r, 'o-', 'LineWidth', 2);
title('(d) Net Return (r)');
xlabel('Tax Rate (%)');
ylabel('Interest Rate');
grid on;

subplot(3, 3, 5);
plot(tax_rates*100, results.TY_ratio, 'o-', 'LineWidth', 2);
title('(e) Transfer to Output Ratio (T/Y)');
xlabel('Tax Rate (%)');
ylabel('T/Y Ratio');
grid on;

subplot(3, 3, 6);
plot(tax_rates*100, results.tax_revenue, 'o-', 'LineWidth', 2);
title('(f) Tax Revenue');
xlabel('Tax Rate (%)');
ylabel('Tax Revenue');
grid on;

subplot(3, 3, 7);
plot(tax_rates*100, results.cons_var, 'o-', 'LineWidth', 2);
title('(g) Consumption Variance');
xlabel('Tax Rate (%)');
ylabel('Variance');
grid on;


path_n_date = [pwd,'\Plots\',datestr(now, 'yyyy-mm-dd-HH-MM'),'\'];
if ~exist(path_n_date, 'dir')
    mkdir(path_n_date);
end
Filename = ['redistribution_analysis_from_',num2str(tax_rates(1)*100),'_to_',num2str(tax_rates(end)*100)];

saveas(gcf, fullfile(path_n_date, [Filename,'.png']));
saveas(gcf, fullfile(path_n_date, [Filename,'.fig']));  % For future edits

% % Save results to a mat file for further analysis
% save('redistribution_results.mat', 'results', 'tax_rates');

% Generate a summary table of results
redistribution_table = table(tax_rates'*100, results.r, results.T, results.Y, results.K, ...
                       results.C, results.TY_ratio, results.tax_revenue, results.cons_var, ...
                       'VariableNames', {'TaxRate_Percent', 'Interest_Rate', 'Transfers', ...
                       'Output', 'Capital', 'Consumption', 'Transfer_Output_Ratio', ...
                       'Tax_Revenue', 'Consumption_Variance'});

% Display the table
disp(redistribution_table);

% Optionally, write to CSV
writetable(redistribution_table, [Filename,'.csv']);
