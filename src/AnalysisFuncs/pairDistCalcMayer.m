% [pDist] = pairDistCalcMayer(n1,n2,n3,l1,l2,lRod)
%
% Calculates the density independent pair correlation function
% g(r,r') = int g(r,r',u,u') du du' = 1 + int fm(r,r',u,u') du du'
%
function [pDist] = pairDistCalcMayer(n1,n2,n3,l1,l2,lRod)
% add paths just in case
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% build phi
dx1 = l1 ./ n1;
dx2 = l2 ./ n2;
dphi = 2 * pi / n3;
x1 = dx1 .* (-n1/2 + 1:n1/2 ); 
x2 = dx2 .* (-n2/2 + 1:n2/2 ); 
phi = 0 : dphi : 2*pi - dphi;
% get MayerFunction
[mayer] = mayerFncHr(n1, n2, n3, l1, l2, lRod) ;
%% intergrate mayer
mayerInt = trapz_periodic( phi, mayer, 3 );
% non-density dependent correlation integrate g(r1,u1,r2,u2)
pDist = 1 + 1 / (2 * pi) * mayerInt;
% center
center1 = round( n1 / 2 ) + 1 ;
center2 = round( n2 / 2 ) + 1 ;
pDist = circshift( circshift( pDist, center1-1, 1 ), center2-1, 2 );
% plot
figure()
imagesc( x1, x2, pDist )
title('pair distribution g(x,y): mayer function averaged')
xlabel('x'); ylabel('y')
colorbar
