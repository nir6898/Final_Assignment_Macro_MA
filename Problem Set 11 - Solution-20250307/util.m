function flowutil = util( c,CRRA_param )
%Utility function - CRRA
%Avoiding -inf in zero
epsilon = 10^(-8);
if CRRA_param==1
    flowutil = log(c);
    flowutil( c < epsilon) = log(epsilon) ;%+ (1/epsilon)*( c(c<epsilon) - epsilon) - (1/(epsilon^2))*0.5*(( c(c<epsilon) - epsilon).^2);
else
    flowutil = (c.^(1-CRRA_param))/(1-CRRA_param);
    flowutil( c < epsilon) = (epsilon.^(1-CRRA_param))/(1-CRRA_param) ;%+ (1/(epsilon^CRRA_param))*( c(c<epsilon) - epsilon ) - (CRRA_param/(epsilon^(CRRA_param + 1)))*0.5*( ( c(c<epsilon) - epsilon ).^2);    
end
end

