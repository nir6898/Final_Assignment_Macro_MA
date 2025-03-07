function val = inv_du( input, CRRA_param )
%The inverse function of the marginal utility of a CRRA
val = 1./(input.^(1/CRRA_param));

%Avoiding -inf in zero
epsilon = 10^(-8);
para = 1/CRRA_param;
val( input < epsilon) = (1./(epsilon.^para)) ;%- para*(1./(epsilon.^(para + 1)))*( input(input<epsilon) - epsilon) + ( para*(para + 1)*(1./(epsilon.^(para + 2))) )*0.5*(( input(input<epsilon) - epsilon).^2);
    
end

