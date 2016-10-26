function [rho] = FixNegDenFnc(rho,systemObj)

% Fix rho if it is negative 

rho_min = min( min ( min( rho ) ) );
epsilon = rho_min / 100;

% keyboard
if rho_min < 0
    rho = rho - (rho_min + epsilon); %Subtract b.c. it's negative
    CurrentNorm = trapz_periodic(x,trapz_periodic(y,trapz_periodic(phi,rho,3),2),1);
    rho = rho .* systemObj.numPart ./ CurrentNorm;
end

end