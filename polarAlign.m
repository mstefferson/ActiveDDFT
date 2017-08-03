% polar aligned potential with no positional dependence
function [v, vFt] = polarAlign( es1,  n3, l3 );
% calculate angles
phi = l3/n3 * (1:n3-1);
% v =  - J u1 dot u2
v = -es1 .* cos( phi )
% Ft
vFt = fftshift( fftn( v ) );
end


