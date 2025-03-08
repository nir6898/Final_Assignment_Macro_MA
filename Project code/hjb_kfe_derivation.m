%% HJB and KFE Equations Derivation
% This file contains the derivation of the HJB and KFE equations for the model
% as required in questions 1-2 of section 2 of the assignment

%% Question 1: HJB Equation

% The HJB equation as a function of (a, zp, IE) is given by:

% For employed households (IE = 1):
% ρV(a, zp, 1) = max_c { u(c) + (∂V/∂a)[ra + zp×w×(1-τ) + T - c] + 
%                         λ_EU[V(a, zp, 0) - V(a, zp, 1)] +
%                         λ_p∑_j[V(a, zp_j, 1) - V(a, zp, 1)]/zp_num }

% For unemployed households (IE = 0):
% ρV(a, zp, 0) = max_c { u(c) + (∂V/∂a)[ra + T - c] + 
%                         λ_UE[V(a, zp, 1) - V(a, zp, 0)] +
%                         λ_p∑_j[V(a, zp_j, 0) - V(a, zp, 0)]/zp_num }

% Where:
% - ρ is the discount rate
% - V(a, zp, IE) is the value function
% - u(c) is the utility function
% - ra + zp×w×(1-τ) + T is the income of employed households
% - ra + T is the income of unemployed households
% - λ_EU is the hazard rate from employment to unemployment
% - λ_UE is the hazard rate from unemployment to employment
% - λ_p is the hazard rate for redrawing the permanent income component
% - zp_num is the number of discretized permanent income states

% The Poisson jump processes for IE and zp explicitly appear in the HJB equation
% through the terms with λ_EU, λ_UE, and λ_p.

%% Question 2: KFE Equation

% The Kolmogorov Forward Equation (KFE) for the distribution g(a, zp, IE, t)
% can be derived as follows:

% Start with the definition of the CDF:
% Prob{a ≤ a', zp = zp', IE = IE'} = ∫_{-∞}^{a'} ∫_{zp = zp'} ∫_{IE = IE'} g(a, zp, IE, t) da dzp dIE

% The evolution of this probability over a small time interval ∆t is given by:
% Prob{a_t+∆t ≤ a', zp_t+∆t = zp', IE_t+∆t = IE'} - Prob{a_t ≤ a', zp_t = zp', IE_t = IE'}

% This change can be decomposed into:
% 1. Drift in assets: ∂g/∂t = -∂/∂a[s(a,zp,IE)×g(a,zp,IE,t)]
% 2. Jumps in employment status:
%    For employed (IE = 1): +λ_UE×g(a,zp,0,t) - λ_EU×g(a,zp,1,t)
%    For unemployed (IE = 0): +λ_EU×g(a,zp,1,t) - λ_UE×g(a,zp,0,t)
% 3. Jumps in permanent income:
%    +λ_p∑_j[g(a,zp_j,IE,t)/zp_num] - λ_p×g(a,zp,IE,t)

% Putting it all together, as ∆t → 0, we get the KFE:

% For employed households (IE = 1):
% ∂g(a,zp,1,t)/∂t = -∂/∂a[s(a,zp,1)×g(a,zp,1,t)] + λ_UE×g(a,zp,0,t) - λ_EU×g(a,zp,1,t)
%                   + λ_p∑_j[g(a,zp_j,1,t)/zp_num] - λ_p×g(a,zp,1,t)

% For unemployed households (IE = 0):
% ∂g(a,zp,0,t)/∂t = -∂/∂a[s(a,zp,0)×g(a,zp,0,t)] + λ_EU×g(a,zp,1,t) - λ_UE×g(a,zp,0,t)
%                   + λ_p∑_j[g(a,zp_j,0,t)/zp_num] - λ_p×g(a,zp,0,t)

% For the stationary distribution, ∂g/∂t = 0, so we have:

% For employed households (IE = 1):
% 0 = -∂/∂a[s(a,zp,1)×g(a,zp,1)] + λ_UE×g(a,zp,0) - λ_EU×g(a,zp,1)
%     + λ_p∑_j[g(a,zp_j,1)/zp_num] - λ_p×g(a,zp,1)

% For unemployed households (IE = 0):
% 0 = -∂/∂a[s(a,zp,0)×g(a,zp,0)] + λ_EU×g(a,zp,1) - λ_UE×g(a,zp,0)
%     + λ_p∑_j[g(a,zp_j,0)/zp_num] - λ_p×g(a,zp,0)

% In matrix form, this can be written as A'g = 0, where A' is the adjoint
% of the operator A that appears in the HJB equation, and g is the vector of
% all distribution masses across the discretized state space.
