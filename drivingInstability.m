rho = linspace( 1.5, 2 );

nem = @nemOrder2;
figure()
plot( rho, nem(rho) )

function s = nemOrder(rho)
  w = 8;
  rhoIN = 1.5;
  s = 1 ./ rho  .* sqrt( 8 / w .* ( rho ./ rhoIN - 1 ) );
end

function s = nemOrder2(rho)
  rhoIN = 1.5;
  f = 2;
  s = 1 ./ (rho) .* f .* ...
    sqrt(  ( f * rho ./ rhoIN - 1 ) .^2 ./ ( rho ./ rhoIN + 1 ) );
end