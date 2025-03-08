% This file contains the utility functions used in the model

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

function inds = z_index_range(zind, knum)
%Z_INDEX_RANGE Finds the indices for a given index of z level
%   Maps from income state index to the corresponding indices in the full state space
    
    firstind = (zind-1)*knum + 1;
    lastind = (zind)*knum;
    inds = firstind:lastind;
end
