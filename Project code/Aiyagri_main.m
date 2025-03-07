clear; clc; close all;
preliminaries;
disp('-------------------- New Run -------------------------')
disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
disp(' ')
% Solver initialization
tic;

%Initializing bisections
net_return_lb =  - params.del;
net_return_ub =  params.rho;
tolK = 10^(-3); % The tolerance
maxit = 40;

for iter = 1:maxit
    net_return = (net_return_lb + net_return_ub)/2;
    [diff, Sol] = model_func( grid, params, net_return);
    if abs(diff) <= tolK
        disp(['Model solved, net returen = ',num2str(net_return)])
        break
    end
    disp(' ')
    disp(['Relative excess demand is = ',num2str(diff),' net return is = ',num2str(net_return)])
    disp(' ')
    %Update the net return
    if diff>0 % excess demand for capital > 0 --> price should be higher
        net_return_lb = net_return;
    else
        net_return_ub = net_return;
    end
end

toc
%% Plot the results:

figure; 
plot(grid.k_vec, Sol.h(1:grid.knum),'DisplayName','High income state','Color','b');hold on; 
plot(grid.k_vec, Sol.h((grid.knum+1):(grid.knum*2)),'DisplayName','Low income state','Color','r'); grid on;
title('Wealth Distribution')
legend

figure
plot(grid.k_vec, Sol.s(1:grid.knum),'DisplayName','High income state','Color','b');hold on; 
plot(grid.k_vec, Sol.s((grid.knum+1):(grid.knum*2)),'DisplayName','Low income state','Color','r'); grid on;
plot(xlim, [0, 0], 'k--', 'DisplayName',"k'-k=0"); % Horizontal line at y = 0
title('Asset Drift')
legend

figure
plot(grid.k_vec, Sol.c(1:grid.knum),'DisplayName','High income state','Color','b');hold on; 
plot(grid.k_vec, Sol.c((grid.knum+1):(grid.knum*2)),'DisplayName','Low income state','Color','r'); grid on;
title('Consumption')
legend


figure
plot(grid.k_vec, Sol.V(1:grid.knum),'DisplayName','High income state','Color','b');hold on; 
plot(grid.k_vec, Sol.V((grid.knum+1):(grid.knum*2)),'DisplayName','Low income state','Color','r'); grid on;
title('Value Function')
legend