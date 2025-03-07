function [diff, Sol] = model_func( grid, params, net_return)

%%=========== Inputs and parameters
rho = params.rho; % Discount rate                    
CRRA = params.CRRA; % Utility function prameter
del = params.del; % Depreciation of capital
alph = params.alph; % Capital share
bet = params.bet; % Labour share

% Grids - unpacking 
znum = params.znum; % Number of human capital levels
zet_vec = params.zet_vec; % Productivity levels

knum = grid.knum;       
k_vec = grid.k_vec;     
dk = grid.dk;  

Zero_fp = 10^(-12); % To determine no drift with finite accuracy



%% ============ Firm behaviour ===========

% Using the analytical firm's problem
gross_return = net_return + del;
K_Dem =  params.Zagg *( alph/gross_return)^( 1/(1-alph) )  ;
w = bet*( (K_Dem/params.Zagg)^alph  );


%%
%================================================
%===============   Households  ==================
%================================================

% Initalize arrays (place holders)
Vn = zeros(znum*knum,1);
V_fd = zeros(size(k_vec));
V_bd = zeros(size(k_vec));

% Drift indicators (place holders)
pos_drift = zeros(size(k_vec));
neg_drift = zeros(size(k_vec));
no_drift = zeros(size(k_vec));

% Return for each level of assets
returns = net_return*k_vec;

% More placeholders
c_upwind = zeros(size(Vn));
u_upwind = zeros(size(Vn));
s_upwind = zeros(size(Vn));

%===== parameters for value function iteration
tol_vf=10^(-8); %Tolerence of value function process

maxit_vf= 5000; %Maximum number of iterations for the value function


flag_vf_converged = 0;
    
VF_Delta_t = 10^6;  

%============Initialization ===================
for zind = 1:znum
    % Initial guess for the relevant value function indices
    Vn(z_index_range(zind, knum)) = util(returns + zet_vec(zind)*w ,CRRA)/(rho); %value from being HtM always
end
%======End initializiation sequence============


% Main loop for household
for vfit = 1:maxit_vf
    %Initializes the LHS matrix - only includes exogenous
    %stochastic process for productivity:
    T_mat = grid.T_mat_base; % -zmat

    %=================== The savings decision ====================
   for zind = 1:znum
        inds = z_index_range(zind, knum);

        %The relevant value function
        V_E = Vn(  inds );

        %Forward differences of the value function        
        V_fd(1:knum-1) = (V_E(2:knum) - V_E(1:knum-1))./dk; 
        V_fd(knum) = du(returns(knum) + zet_vec(zind)*w,CRRA); % edge case

        %Backwards differences of the value function
        V_bd(2:knum) = (V_E(2:knum) - V_E(1:knum-1))./dk; 
        V_bd(1) = du(returns(1) + zet_vec(zind)*w ,CRRA); % edge case

        %Policy functions using forward differences of the value function
        c_fd = inv_du(V_fd, CRRA);

        %Policy functions using backwards differences of the value function
        c_bd = inv_du(V_bd, CRRA);

        %Saving functions using forward differences of the value function
        s_fd = returns + zet_vec(zind)*w - c_fd ; 

        %Saving functions using backwards differences of the value function
        s_bd = returns + zet_vec(zind)*w - c_bd; 

        %If savings are zero:
        c_nd = returns + zet_vec(zind)*w ;
        V_nd=du(c_nd,CRRA);

        %indicators for drifts:        
        pos_drift = s_fd>Zero_fp; % Zero_fp ~ 0. it's used to avoid floating point imprecisions
        neg_drift = s_bd<-Zero_fp;
        no_drift = 1 - pos_drift - neg_drift ;

        %derivatives
        dV_upwind = pos_drift.*V_fd + neg_drift.*V_bd + no_drift.*V_nd;

        c_upwind(  inds  ) = inv_du(dV_upwind,CRRA);
        u_upwind(  inds  ) = util(c_upwind(  inds  ),CRRA);
        s_upwind(  inds  ) = returns + zet_vec(zind)*w - c_upwind(  inds  );


        %Constructing the D matrix in the scheme for each value
        %function. Nodes are ordered by the indices of the k_vec,
        %lower to higher. This matrix represents the upwind opperator
        %using the decision rules for assets - forward difference for
        %positive state drift and backward difference for
        %negtive drift (backwards in time)

        main_diag = -[max(s_fd(1:(knum-1)),0)./dk;0] + [0;min(s_bd(2:knum),0)./dk]; 
        lower_diag = -min(s_bd(2:knum),0)./dk;
        upper_diag = max(s_fd(1:knum-1),0)./dk;
        D = spdiags(main_diag,0,knum,knum) + spdiags(lower_diag,-1,knum,knum) + spdiags([0;upper_diag],1,knum,knum);
        T_mat(  inds ,  inds  ) = T_mat(  inds ,  inds  ) -D;

   end
    

    AAA = grid.T_mat_III .*((1/VF_Delta_t) + rho); % The last part of the matrix we are building
    T_mat = T_mat + AAA;


    %========== The full discretization scheme:
    %RHS of the numerical scheme
    t_vector = u_upwind + Vn/VF_Delta_t;
    %Value function
    Vnp1 = T_mat\t_vector;		
    %Stopping criterion
    dist_V = max(abs(Vnp1-Vn));		
    
    %Prepare for the next iteration 
    Vn = Vnp1;

    %convergence check
    if dist_V<tol_vf
        disp(['VF solved. Error was =',num2str(dist_V)])
        flag_vf_converged = 1;
        break
    end        
    
end
        
   

%%
%==============================================================
%=======Distributions==========================================
%==============================================================

%Between the HACT algorithm and my implementation - T_mat is the adjoint
%transpose operator but on the LHS rather than the RHS thus I transpose and
%flip the sign.

% We built the martix M earlier, which has the matrix A as a component.
% We want A back:
AAA = grid.T_mat_III .*((1/VF_Delta_t) + rho);
A = - (T_mat - AAA);

% The matrix that governs the distribution is A transposed:
G = transpose(A);
Sol.KFE_mat = G;

%The infinitly lived version with a known pop mass, the KFE is solved by:
% Stationary distribution, all mass changes are 0:
c = zeros(znum*knum, 1);

% But in order to have a unique solution, we change the first "equation" to be: the sum of masses=1
c(1) = 1;
G(1,:) = 1;

%distribution:

h = G\c;
% sum(h) should be 1, and sum(Sol.h(1:3001)) should be 0.5

%==================Storing policies===================
Sol.V = Vn;
Sol.c = c_upwind;
Sol.s = s_upwind;
Sol.h = h;

%================Aggregation========================

%Aggregate assets

% The aggregate supply of assets:
K_sup = sum(repmat(k_vec,znum,1).*h); %rep mat(vec, n, 1) just concatinates the vector to itself n times 
Sol.Kagg = K_sup;

% Aggregate consumption:
Sol.Agg_Cons = sum(c_upwind.*h);

% Aggregate output:
Sol.Agg_output = (K_sup^alph)*(params.Zagg^bet);

%% Finishing:

Sol.K_Dem = K_Dem;

%Captial market clearing condition:
diff = (K_Dem - K_sup)/K_sup;

end

