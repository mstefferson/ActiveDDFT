function [rho] = FixNegDenFnc(rho)

% Fix rho if it is negative 

rho_min = min( min ( min( rho ) ) );
epsilon = rho_min / 100;

% keyboard
if rho_min < 0
    
    rho = rho - (rho_min + epsilon); %Subtract b.c. it's negative
    
end

end