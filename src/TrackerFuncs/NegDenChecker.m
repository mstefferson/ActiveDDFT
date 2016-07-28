function NegDenChecker(rho)

if min(min(min(rho))) < 0
   error('Density is negative') 
end


end
