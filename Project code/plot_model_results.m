function plot_model_results(grid_struct, params, Sol)
% PLOT_MODEL_RESULTS Generates all the required plots for the model results
%   This function creates the plots required in section 2.5 of the assignment

path_n_date = [pwd,'\Plots\',datestr(now, 'yyyy-mm-dd-HH-MM'),'\'];
if ~exist(path_n_date, 'dir')
    mkdir(path_n_date);
end

% Extract grid and parameters for convenience
k_vec = grid_struct.k_vec;
knum = grid_struct.knum;
zp_num = params.zp_num;

% Select specific zp grid points to plot as requested
plot_zp_indices = [1, 3, 7, 10]; % As required in the assignment

%% Figure 1: Consumption policy functions of employed households
figure('Position', [100, 100, 800, 600]);
colors = jet(length(plot_zp_indices));

for i = 1:length(plot_zp_indices)
    zp_idx = plot_zp_indices(i);
    z_idx = zp_idx; % Employed state
    inds = ((z_idx-1)*knum + 1):((z_idx)*knum);
    
    plot(k_vec, Sol.c(inds), 'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', ['z_p = ', num2str(params.zp_grid(zp_idx), '%.2f')]);
    hold on;
end

title('Consumption Policy Functions - Employed Households');
xlabel('Assets (a)');
ylabel('Consumption (c)');
grid on;
legend('Location', 'Southeast');
saveas(gcf, fullfile(path_n_date, 'employed_consumption.png'));

%% Figure 2: Consumption policy functions of unemployed households
figure('Position', [100, 100, 800, 600]);

for i = 1:length(plot_zp_indices)
    zp_idx = plot_zp_indices(i);
    z_idx = zp_idx + zp_num; % Unemployed state with same zp
    inds = ((z_idx-1)*knum + 1):((z_idx)*knum);
    
    plot(k_vec, Sol.c(inds), 'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', ['z_p = ', num2str(params.zp_grid(zp_idx), '%.2f')]);
    hold on;
end

title('Consumption Policy Functions - Unemployed Households');
xlabel('Assets (a)');
ylabel('Consumption (c)');
grid on;
legend('Location', 'Southeast');
saveas(gcf, fullfile(path_n_date, 'unemployed_consumption.png'));

%% Figure 3: Sparsity structure of the A matrix
figure('Position', [100, 100, 800, 800]);

% Get the A matrix from the final iteration
AAA = grid_struct.T_mat_III .* ((1/params.VF_Delta_t) + params.rho);
A = -(Sol.KFE_mat' - AAA);

% Plot the sparsity pattern
spy(A);
title('Sparsity Structure of the A Matrix from HJB Discretization');
xlabel('Column Index');
ylabel('Row Index');
saveas(gcf, fullfile(path_n_date, 'A_matrix_sparsity.png'));

%% Figure 4: Stationary distribution of assets
figure('Position', [100, 100, 800, 600]);

% Aggregate distribution over all states
agg_dist = zeros(knum, 1);
for i = 1:params.znum
    inds = ((i-1)*knum + 1):((i)*knum);
    agg_dist = agg_dist + Sol.h(inds);
end

% Plot as histogram
bar(k_vec, agg_dist, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none');
title('Stationary Distribution of Assets in the Population');
xlabel('Assets (a)');
ylabel('Density');
grid on;
xlim([0, min(grid_struct.k_max, 50)]); % Limit x-axis for better visibility
saveas(gcf, fullfile(path_n_date, 'asset_distribution.png'));

%% Figure 5: Net contribution by zp level
figure('Position', [100, 100, 800, 600]);

net_contribution = zeros(params.zp_num, 1);
zp_values = params.zp_grid;

% Calculate net contribution for each zp level
for i = 1:params.zp_num
    % Employed with this zp
    e_inds = ((i-1)*knum + 1):((i)*knum);
    e_mass = sum(Sol.h(e_inds));
    e_contrib = params.tau * params.zp_grid(i) * Sol.wage * e_mass;
    
    % Unemployed with this zp
    u_inds = ((i + params.zp_num-1)*knum + 1):((i + params.zp_num)*knum);
    u_mass = sum(Sol.h(u_inds));
    
    % Total mass with this zp
    total_mass = e_mass + u_mass;
    
    % Net contribution (taxes paid minus transfers received)
    net_contribution(i) = e_contrib - Sol.T * total_mass;
end

% Create bar chart
bar(1:params.zp_num, net_contribution, 'FaceColor', [0.4 0.7 0.3]);
title('Net Contribution by Income Level (τwz - T)');
xlabel('Permanent Income Level (z_p)');
ylabel('Net Contribution');
xticks(1:params.zp_num);
xticklabels(arrayfun(@(x) sprintf('%.2f', x), zp_values, 'UniformOutput', false));
xtickangle(45);
grid on;

% Add a horizontal line at zero for reference
hold on;
plot(xlim, [0, 0], 'r--', 'LineWidth', 1.5);
legend('Net Contribution', 'Zero Line', 'Location', 'Southeast');
saveas(gcf, fullfile(path_n_date, 'net_contribution.png'));

%% Additional useful plots

% Savings policy function (asset drift)
figure('Position', [100, 100, 800, 600]);
for i = 1:length(plot_zp_indices)
    % Employed
    zp_idx = plot_zp_indices(i);
    z_idx = zp_idx;
    inds = ((z_idx-1)*knum + 1):((z_idx)*knum);
    
    plot(k_vec, Sol.s(inds), 'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', ['Employed, z_p = ', num2str(params.zp_grid(zp_idx), '%.2f')]);
    hold on;
    
    % Unemployed
    z_idx = zp_idx + zp_num;
    inds = ((z_idx-1)*knum + 1):((z_idx)*knum);
    
    plot(k_vec, Sol.s(inds), 'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', '--', ...
         'DisplayName', ['Unemployed, z_p = ', num2str(params.zp_grid(zp_idx), '%.2f')]);
end

% Add a horizontal line at zero drift
plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero Drift');

title('Asset Drift (s = da/dt)');
xlabel('Assets (a)');
ylabel('Savings (da/dt)');
grid on;
legend('Location', 'Best');
saveas(gcf, fullfile(path_n_date, 'asset_drift.png'));

% Value function
figure('Position', [100, 100, 800, 600]);
for i = 1:length(plot_zp_indices)
    % Employed
    zp_idx = plot_zp_indices(i);
    z_idx = zp_idx;
    inds = ((z_idx-1)*knum + 1):((z_idx)*knum);
    
    plot(k_vec, Sol.V(inds), 'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', ['Employed, z_p = ', num2str(params.zp_grid(zp_idx), '%.2f')]);
    hold on;
    
    % Unemployed
    z_idx = zp_idx + zp_num;
    inds = ((z_idx-1)*knum + 1):((z_idx)*knum);
    
    plot(k_vec, Sol.V(inds), 'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', '--', ...
         'DisplayName', ['Unemployed, z_p = ', num2str(params.zp_grid(zp_idx), '%.2f')]);
end

title('Value Functions');
xlabel('Assets (a)');
ylabel('Value');
grid on;
legend('Location', 'Southeast');
saveas(gcf, fullfile(path_n_date, 'value_functions.png'));
end