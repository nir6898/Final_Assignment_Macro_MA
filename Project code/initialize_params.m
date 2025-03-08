%% Initialize parameters and grids for the model
% This script creates two structs:
% 1. params: Contains all model parameters
% 2. grid: Contains grids and pre-initialized matrices

%% Household parameters
params.rho = 0.1;             % Discount rate
params.borrowing_lim = 0;     % Borrowing limit (a = 0)
params.del = 0.048;           % Depreciation of capital
params.CRRA = 2;              % Utility function parameter
params.alph = 0.36;           % Capital share
params.bet = 1 - params.alph; % Labour share
params.tau = 0.4;             % Labor income tax rate

%% Income process parameters
% Define employment/unemployment process
params.lambda_EU = 0.4;       % Employed to unemployed transition rate
params.lambda_UE = 5.4;       % Unemployed to employed transition rate

% Expected duration calculations
expected_employment_duration = 1/params.lambda_EU;
expected_unemployment_duration = 1/params.lambda_UE;
fprintf('Expected employment duration: %.2f years\n', expected_employment_duration);
fprintf('Expected unemployment duration: %.2f years\n', expected_unemployment_duration);

% Define permanent income component
params.zeta = 1.5;            % Pareto distribution parameter
params.zp_num = 10;           % Number of grid points for permanent income

% Calibrate permanent income shock rate
params.lambda_p = 1/30;       % Expected duration of 30 years

% Create discretized Pareto distribution
min_zp = 1;  % Minimum value for zp
zp_grid = zeros(params.zp_num, 1);

% Using quantile function to discretize Pareto distribution
for i = 1:params.zp_num
    p = (i - 0.5) / params.zp_num;  % Midpoint of probability mass
    zp_grid(i) = min_zp * (1 - p)^(-1/params.zeta);
end

% Normalize zp grid so the mean is 1
zp_mean = mean(zp_grid);
zp_grid = zp_grid / zp_mean;

params.zp_grid = zp_grid;

% Create the full income state space (zp × IE)
params.znum = 2 * params.zp_num;  % Total number of income states
params.zet_vec = zeros(params.znum, 1);

% Fill the income vector: first half is employed states, second half is unemployed
for i = 1:params.zp_num
    params.zet_vec(i) = params.zp_grid(i);                 % Employed with zp(i)
    params.zet_vec(i + params.zp_num) = 0;                 % Unemployed with zp(i)
end

% Create transition matrix for the income process
M_Z = zeros(params.znum, params.znum);

% Fill transition matrix
for i = 1:params.zp_num
    for j = 1:params.zp_num
        % From employed with zp(i) to employed with zp(j)
        if i == j  % No change in zp
            M_Z(i, j) = -params.lambda_EU - params.lambda_p;
        else       % Change in zp (redrawn from distribution)
            M_Z(i, j) = params.lambda_p * (1/params.zp_num);
        end
        
        % From employed with zp(i) to unemployed with zp(j)
        if i == j  % Same zp, but becomes unemployed
            M_Z(i, j + params.zp_num) = params.lambda_EU;
        else       % Different zp (this shouldn't happen in this model)
            M_Z(i, j + params.zp_num) = 0;
        end
        
        % From unemployed with zp(i) to employed with zp(j)
        if i == j  % Same zp, but becomes employed
            M_Z(i + params.zp_num, j) = params.lambda_UE;
        else       % Different zp (this shouldn't happen in this model)
            M_Z(i + params.zp_num, j) = 0;
        end
        
        % From unemployed with zp(i) to unemployed with zp(j)
        if i == j  % No change in zp
            M_Z(i + params.zp_num, j + params.zp_num) = -params.lambda_UE - params.lambda_p;
        else       % Change in zp (redrawn from distribution)
            M_Z(i + params.zp_num, j + params.zp_num) = params.lambda_p * (1/params.zp_num);
        end
    end
end

params.M_Z = M_Z;

% Calculate stationary distribution of income states
% Using the Kolmogorov Forward Equation for the Poisson jump process
KF = params.M_Z';
% Add the condition that probabilities sum to 1
KF(1,:) = 1;
z_dist = KF\[1; zeros(params.znum-1, 1)];

% Normalize to ensure the sum is exactly 1 (numerical precision)
z_dist = z_dist / sum(z_dist);

% Calculate aggregate labor supply for initial wage determination
params.Zagg = sum(z_dist .* params.zet_vec);
fprintf('Aggregate labor supply: %.4f\n', params.Zagg);

% Store the stationary distribution
params.z_dist = z_dist;
params.Agg_Labor = params.Zagg;

%% Asset grid construction
grid.k_max = 300;                      % Maximum asset level
grid.k_min = -params.borrowing_lim;    % Minimum asset level (no borrowing)
grid.knum = 1000;                      % Number of grid points
grid.k_vec = linspace(grid.k_min, grid.k_max, grid.knum)'; % Asset grid
grid.dk = grid.k_vec(2) - grid.k_vec(1);  % Grid step size

%% Construct block matrices for HJB and KFE equations
% Create sparse identity matrices
grid.III = speye(grid.knum, grid.knum);
state_num = params.znum * grid.knum;
grid.T_mat_III = speye(state_num, state_num);

% Create transition matrix for income process
zmat = sparse(state_num, state_num);
trans_mat = params.M_Z;

% Populate matrix representing income transitions
for zind = 1:params.znum
    for zprim = 1:params.znum
        block_inds_z = ((zind-1)*grid.knum + 1):((zind)*grid.knum);
        block_inds_zp = ((zprim-1)*grid.knum + 1):((zprim)*grid.knum);
        
        if trans_mat(zind, zprim) ~= 0
            zmat(block_inds_z, block_inds_zp) = trans_mat(zind, zprim) * grid.III;
        end
    end
end

% Base transition matrix (only income transitions, no asset drift yet)
grid.T_mat_base = -zmat;

% Parameters for numerical solution
params.tol_vf = 1e-8;      % Value function convergence tolerance
params.maxit_vf = 5000;    % Maximum value function iterations
params.VF_Delta_t = 1e6;   % Time step for implicit scheme
params.Zero_fp = 1e-12;    % Tolerance for zero drift
