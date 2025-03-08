function [diff, Sol] = model_func(grid, params, net_return, T)
% MODEL_FUNC Implements steps 2-5 of the algorithm
%   Solves the firm's problem, household's problem, and finds the stationary
%   distribution for given interest rate and transfer
%
%   Inputs:
%     grid:       Struct with grids and matrices
%     params:     Struct with model parameters
%     net_return: Net return on capital (r)
%     T:          Lump-sum transfer
%
%   Outputs:
%     diff:       Capital market clearing error
%     Sol:        Struct with solution objects (policies, distribution, etc.)

%%=========== Inputs and parameters
rho = params.rho;              % Discount rate                    
CRRA = params.CRRA;            % Utility function parameter
del = params.del;              % Depreciation of capital
alph = params.alph;            % Capital share
bet = params.bet;              % Labour share
tau = params.tau;              % Labor income tax rate

% Grids - unpacking 
znum = params.znum;            % Number of income states
zet_vec = params.zet_vec;      % Income levels

knum = grid.knum;       
k_vec = grid.k_vec;     
dk = grid.dk;  

Zero_fp = params.Zero_fp;      % To determine no drift with finite accuracy

%% ============ Firm behaviour ===========

% Using the analytical firm's problem
gross_return = net_return + del;
K_Dem = params.Zagg * (alph/gross_return)^(1/(1-alph));
w = bet * ((K_Dem/params.Zagg)^alph);

%%
%================================================
%===============   Households  ==================
%================================================

% Initialize arrays (place holders)
Vn = zeros(znum*knum, 1);
V_fd = zeros(size(k_vec));
V_bd = zeros(size(k_vec));

% Drift indicators (place holders)
pos_drift = zeros(size(k_vec));
neg_drift = zeros(size(k_vec));
no_drift = zeros(size(k_vec));

% Return for each level of assets
returns = net_return*k_vec;

% More placeholders
c_upwind = zeros(znum*knum, 1);
u_upwind = zeros(znum*knum, 1);
s_upwind = zeros(znum*knum, 1);

%===== parameters for value function iteration
tol_vf = params.tol_vf;        % Tolerance of value function process
maxit_vf = params.maxit_vf;    % Maximum number of iterations for the value function
VF_Delta_t = params.VF_Delta_t;  

flag_vf_converged = 0;

%============Initialization ===================
for zind = 1:znum
    % Initial guess for the relevant value function indices
    inds = ((zind-1)*knum + 1):((zind)*knum);
    
    % Calculate after-tax income for the initial guess
    if zind <= params.zp_num % Employed states
        income = zet_vec(zind) * w * (1-tau) + T;
    else                     % Unemployed states
        income = T;
    end
    
    % Initial guess: value from consuming all income forever
    Vn(inds) = util(returns + income, CRRA)/(rho);
end
%======End initialization sequence============

% Main loop for household
for vfit = 1:maxit_vf
    % Initializes the LHS matrix - only includes exogenous
    % stochastic process for productivity:
    T_mat = grid.T_mat_base;  % -zmat

    %=================== The savings decision ====================
    for zind = 1:znum
        inds = ((zind-1)*knum + 1):((zind)*knum);

        % Calculate after-tax income
        if zind <= params.zp_num  % Employed states
            income = zet_vec(zind) * w * (1-tau) + T;
        else                      % Unemployed states
            income = T;
        end

        % The relevant value function
        V_E = Vn(inds);

        % Forward differences of the value function        
        V_fd(1:knum-1) = (V_E(2:knum) - V_E(1:knum-1))./dk; 
        V_fd(knum) = du(returns(knum) + income, CRRA);  % edge case

        % Backwards differences of the value function
        V_bd(2:knum) = (V_E(2:knum) - V_E(1:knum-1))./dk; 
        V_bd(1) = du(returns(1) + income, CRRA);  % edge case

        % Policy functions using forward differences of the value function
        c_fd = inv_du(V_fd, CRRA);

        % Policy functions using backwards differences of the value function
        c_bd = inv_du(V_bd, CRRA);

        % Saving functions using forward differences of the value function
        s_fd = returns + income - c_fd; 

        % Saving functions using backwards differences of the value function
        s_bd = returns + income - c_bd; 

        % If savings are zero:
        c_nd = returns + income;
        V_nd = du(c_nd, CRRA);

        % Indicators for drifts:        
        pos_drift = s_fd > Zero_fp;  % Zero_fp ~ 0, used to avoid floating point imprecisions
        neg_drift = s_bd < -Zero_fp;
        no_drift = 1 - pos_drift - neg_drift;

        % Derivatives
        dV_upwind = pos_drift.*V_fd + neg_drift.*V_bd + no_drift.*V_nd;

        c_upwind(inds) = inv_du(dV_upwind, CRRA);
        u_upwind(inds) = util(c_upwind(inds), CRRA);
        s_upwind(inds) = returns + income - c_upwind(inds);

        % Constructing the D matrix in the scheme for each value
        % function. Nodes are ordered by the indices of the k_vec,
        % lower to higher. This matrix represents the upwind operator
        % using the decision rules for assets - forward difference for
        % positive state drift and backward difference for
        % negative drift (backwards in time)
        main_diag = -[max(s_fd(1:(knum-1)),0)./dk;0] + [0;min(s_bd(2:knum),0)./dk]; 
        lower_diag = -min(s_bd(2:knum),0)./dk;
        upper_diag = max(s_fd(1:knum-1),0)./dk;
        D = spdiags(main_diag,0,knum,knum) + spdiags(lower_diag,-1,knum,knum) + spdiags([0;upper_diag],1,knum,knum);
        T_mat(inds, inds) = T_mat(inds, inds) - D;
    end

    AAA = grid.T_mat_III .* ((1/VF_Delta_t) + rho);  % The last part of the matrix we are building
    T_mat = T_mat + AAA;

    % The full discretization scheme:
    % RHS of the numerical scheme
    t_vector = u_upwind + Vn/VF_Delta_t;
    % Value function
    Vnp1 = T_mat\t_vector;		
    % Stopping criterion
    dist_V = max(abs(Vnp1-Vn));		
    
    % Prepare for the next iteration 
    Vn = Vnp1;

    % Convergence check
    if dist_V < tol_vf
        % disp(['VF solved. Error was = ', num2str(dist_V)]);
        flag_vf_converged = 1;
        break;
    end
end

%%
%==============================================================
%=======Distributions==========================================
%==============================================================

% Between the HACT algorithm and this implementation - T_mat is the adjoint
% transpose operator but on the LHS rather than the RHS thus we transpose and
% flip the sign.

% We built the matrix M earlier, which has the matrix A as a component.
% We want A back:
AAA = grid.T_mat_III .* ((1/VF_Delta_t) + rho);
A = -(T_mat - AAA);

% The matrix that governs the distribution is A transposed:
G = A';
Sol.KFE_mat = G;

% The infinitely lived version with a known pop mass, the KFE is solved by:
% Stationary distribution, all mass changes are 0:
c = zeros(znum*knum, 1);

% But in order to have a unique solution, we change the first "equation" to be: the sum of masses=1
c(1) = 1;
G(1,:) = 1;

% Distribution:
h = G\c;

% Calculate total employed and unemployed mass
employed_mass = sum(h(1:(params.zp_num*knum)));
unemployed_mass = sum(h((params.zp_num*knum+1):end));
% fprintf('Employed mass: %.4f, Unemployed mass: %.4f\n', employed_mass, unemployed_mass);

% Check that the distribution sums to approximately 1
% disp(['Sum of distribution: ', num2str(sum(h))]);

%==================Storing policies===================
Sol.V = Vn;
Sol.c = c_upwind;
Sol.s = s_upwind;
Sol.h = h;

%================Aggregation========================

% Aggregate assets (capital supply)
K_sup = 0;
for zind = 1:znum
    inds = ((zind-1)*knum + 1):((zind)*knum);
    K_sup = K_sup + sum(k_vec .* h(inds));
end
Sol.Kagg = K_sup;

% Aggregate consumption
Sol.Agg_Cons = sum(c_upwind .* h);

% Aggregate labor (employed workers weighted by productivity)
Sol.Agg_Labor = params.Zagg;

% Calculate aggregate labor income (taxable income)
Sol.Agg_Labor_Income = 0;
for zind = 1:params.zp_num  % Only employed states
    inds = ((zind-1)*knum + 1):((zind)*knum);
    Sol.Agg_Labor_Income = Sol.Agg_Labor_Income + sum(zet_vec(zind) * h(inds));
end
Sol.Agg_Labor_Income = Sol.Agg_Labor_Income * w;

% Aggregate output using production function
Sol.Agg_output = (K_sup^alph) * (params.Zagg^bet);

% Store other useful information
Sol.K_Dem = K_Dem;
Sol.wage = w;
Sol.net_return = net_return;
Sol.T = T;

% Capital market clearing condition (error)
diff = (K_Dem - K_sup)/K_sup;

end

% Utility functions
function flowutil = util(c, CRRA_param)
%UTIL Utility function - CRRA
%   Computes the utility for consumption c with CRRA parameter
    
    % Avoid -inf in zero by setting a small epsilon
    epsilon = 1e-8;
    
    if CRRA_param == 1
        % Log utility case
        flowutil = log(c);
        flowutil(c < epsilon) = log(epsilon);
    else
        % General CRRA case
        flowutil = (c.^(1-CRRA_param))/(1-CRRA_param);
        flowutil(c < epsilon) = (epsilon.^(1-CRRA_param))/(1-CRRA_param);
    end
end

function val = du(input, CRRA_param)
%DU Marginal utility of a CRRA utility function
%   Computes the first derivative of the utility function
    
    % Standard CRRA marginal utility
    val = input.^(-CRRA_param);
    
    % Avoid -inf in zero
    epsilon = 1e-8;
    val(input < epsilon) = epsilon^(-CRRA_param);
end

function val = inv_du(input, CRRA_param)
%INV_DU Inverse function of the marginal utility of a CRRA
%   Computes the consumption that yields a given marginal utility
    
    % Standard inverse of CRRA marginal utility
    val = input.^(-1/CRRA_param);
    
    % Avoid issues near zero
    epsilon = 1e-8;
    para = 1/CRRA_param;
    val(input < epsilon) = epsilon^(-para);
end