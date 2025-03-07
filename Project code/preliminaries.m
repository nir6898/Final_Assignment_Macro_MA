%Household parameters
params.rho = (1-0.95)/0.95;
params.borrowing_lim = 0; 
params.del = 0.05; % Depreciation of capital
params.CRRA = 1; %Utility function prameter
params.alph = 0.34; %Capital share
params.bet =  0.66; %Labour share - it's just 1 - alpha


%===============Labour income process
% ==== Stochastic process for labour income ==============
params.zet_vec = [3 ; 1];
params.znum = length(params.zet_vec);

params.M_Z = [-0.4 0.4; 0.4 -0.4];

%Using the KFE for the Poisson jump process to solve for the masses of each
%household's type
KF = params.M_Z';
KF(1,:) = 1; % Imposing sum of masses to be 1
z_dist = KF\[1;0]; % Finding the stationary aggregate wealth distribution

params.Zagg = sum(z_dist.*params.zet_vec);

%% Grids construction
%========Grid - asset=======

grid.k_max = 60 ;
grid.k_min = -params.borrowing_lim ;
grid.knum = 3001;
k_range = grid.k_max - grid.k_min;
k_vec = linspace(-params.borrowing_lim, grid.k_max, grid.knum)';
dk = k_vec(2:grid.knum) - k_vec(1:(grid.knum-1));
grid.k_vec = k_vec;
grid.dk = dk;  



%========== Construct block matrices for HJB block ========
znum = params.znum;
grid.III = speye(grid.knum,grid.knum); % Sparse identity matrix
state_num = (znum)*grid.knum; 

%Preallocating memory for the the transition matrix
T_mat_III = speye(state_num,state_num); % Sparse identity matrix
grid.T_mat_III = T_mat_III;

%Inflate the labour income process to the size of the asset grid 

zmat = speye(state_num,state_num);
trans_mat = params.M_Z;
% zmat represents the infinitesimal generator in our linear system
for zind = 1:znum
    for zprim = 1:znum
        zmat( ((zind-1)*grid.knum + 1):((zind)*grid.knum), ((zprim-1)*grid.knum + 1):((zprim)*grid.knum) ) =  trans_mat(zind,zprim)*grid.III;
    end
end




%% Preinitializing the transition matrix (exogenous parts only) ========
%=======================The stochastic process for z ============================= 

% constructing T_mat as if it is already on the LHS of the scheme

T_mat_base = - zmat;

grid.T_mat_base = T_mat_base;

%save('gridk.mat','grid')
%save('calib.mat','params')