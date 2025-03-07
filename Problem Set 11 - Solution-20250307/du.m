
function  val = du(  input, CRRA_param  )
%The marginal utility of a CRRA
val = 1./(input.^CRRA_param);
%Avoiding -inf in zero
epsilon = 10^(-8);
val( input < epsilon) = (1./(epsilon.^CRRA_param)) ;%- CRRA_param*(1./(epsilon.^(CRRA_param + 1)))*( input(input<epsilon) - epsilon) + (CRRA_param*(CRRA_param + 1)*(1./(epsilon.^(CRRA_param + 2))) )*0.5*(( input(input<epsilon) - epsilon).^2);
end

