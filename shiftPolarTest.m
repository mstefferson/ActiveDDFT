% 1D example

Nm  = 32;
dphi = 2 * pi / Nm;

phi = 0 : dphi : 2*pi - dphi;

% 1D distro
f = 1 + cos( 2 .* phi ); 
f(1) = f(1) + 0.01;
normTemp = trapz_periodic( phi, f );
f =  1 ./ normTemp .* f;

% perturbation
pAmp = 1 / (2*pi);
pVar = 2 * pi / 2;

% p = aphi .* cos(phi);
[~, maxfInd] = max( f );

p =  gaussShift( phi, pAmp, pVar, phi( maxfInd ) );

% New distro
gNoSpat = p + f;
gmin = min(gNoSpat);
if gmin < 0
  gNoSpat = gNoSpat - gmin;
end
normTemp = trapz_periodic( phi, gNoSpat );
gNoSpat =  1 ./ normTemp .* gNoSpat;

figure()
plot( phi, f, phi, p, phi, gNoSpat )

legend( 'original', 'perturbation', 'final')

%% 2D example

% grid stuff
c   = 2;
Nx  = 24;
Lx  = 1;
Nm  = 24;
dphi = 2 * pi / Nm;
dx = Lx / Nx;

% perturbation
pAmpPhi = 1 / (2*pi); % amptiude of gaussian for phi as function of x
pVarPhi = Lx / 2;  % variance of amps in pos
pVar = 2 * pi / 4;
phi = 0 : dphi : 2*pi - dphi;
x = 0 : dx : Lx - dx;

% [phi2, x2] = meshgrid( phi,x );

% 1D distro
f = 1 + cos( 2 .* phi );
normTemp = trapz_periodic( phi, f );
f =  1 ./ normTemp .* f;
f = repmat( f, [Nx,1] );

rho = c .* f;

% perturbation


% ampitude for gaussian that changes vs position
aX  =  gaussShift( x, pAmpPhi, Lx / 2, pVarPhi ); % amp(x)
% p2 = zeros(Nx,Nm);

[~, maxfInd] = max( rho(1,:) );
p = gaussShift( phi, 1, pVar, phi( maxfInd ) );

p2 = repmat( p, [Nx, 1] );

aX2 =  repmat( aX', [1 Nm] );

gNoSpat =rho + aX2 .* p2;

% 
% % homogenous
% for ii = 1:Nx
%   [~, maxfInd] = max( rho(ii,:) );
%   p2( ii, : ) = gaussShift( phi, aX(ii), pVar, phi( maxfInd ) );
%   gNoSpat( ii, : ) = p2( ii, : ) + rho( ii,: );
%   
%   normTemp = trapz_periodic( phi, gNoSpat( ii, : ) );
%   gNoSpat( ii, : ) = normTemp ./ c .* gNoSpat( ii, : ) ;
% end

gNoSpat = p2 + rho;
gmin = min( min(gNoSpat) );

if gmin < 0
  gNoSpat = gNoSpat - gmin;
end

normTemp = trapz_periodic( x, trapz_periodic( phi, gNoSpat, 2 ), 1 );
gNoSpat =  c  ./ normTemp .* gNoSpat;

figure()
subplot(1,3,1)
imagesc( x, phi, rho')
xlabel('x'); ylabel( 'phi');
title( 'orig' )
colorbar
axis square

subplot(1,3,2)
imagesc( x, phi, gNoSpat')
xlabel('x'); ylabel( 'phi');
title('no conc variance')
colorbar
axis square

% subplot(2,2,3)
% imagesc( x, phi, ax')
% xlabel('x'); ylabel( 'phi');
% title('spatial slope')
% colorbar
% axis square
% 
% subplot(2,2,4)
% imagesc( x, phi, p')
% xlabel('x'); ylabel( 'phi');
% title('pert')
% colorbar
% axis square

legend( 'original', 'perturbation', 'final')

%% 3D
normTot = 2;
Nx  = 24;
Lx  = 1;
Nm  = 24;
dphi = 2 * pi / Nm;
dx = Lx / Nx;

phi = 0 : dphi : 2*pi - dphi;
x = 0 : dx : Lx - dx;

[phi2, x2] = meshgrid( phi,x );

% 1D distro
f = 1 + cos( 2 .* phi );
normTemp = trapz_periodic( phi, f );
f =  1 ./ normTemp .* f;
f = repmat( f, [Nx,1] );

% perturbation
aphi = 1 / (2*pi);
aspat = 1/ (2*pi);
ax = -aphi .* x2  .* (x2 -Lx) ;

pnospat = ax .* cos(phi2);
gVar = Lx/4;
pSpat = ax .* cos(phi2) + aspat .* exp( - ( x2 - Lx/2 ) .^ 2  ./ ( 2 * gVar .^2 )  );
pNoSpat = ax .* cos(phi2);

gNoSpat = f + pNoSpat ;
gSpat = f + pSpat;

% Fix negatives
gmin = min( min(gNoSpat) );
if gmin < 0
  gNoSpat = gNoSpat - gmin;
end

gmin = min( min(gSpat) );
if gmin < 0
  gSpat = gSpat - gmin;
end

% g inhomogenous (in x)
normTemp = trapz_periodic( trapz_periodic( phi, gSpat, 2 ), 1 );
ginhom =  normTot  ./ normTemp .* gSpat;
normTemp = trapz_periodic( phi, ginhom, 2 );

% g homogenous (in x)
ghom = gSpat;
normTemp = trapz_periodic( phi, gSpat, 2 );

for ii = 1:Nx
  ghom(ii,:) =  normTot  ./ normTemp(ii) .* gSpat(ii,:);
end

normTemp = trapz_periodic( phi, ghom, 2 );

figure()
subplot(2,2,1)
imagesc( x, phi, f')
xlabel('x'); ylabel( 'phi');
title( 'orig' )
colorbar
axis square

subplot(2,2,2)
imagesc( x, phi, gNoSpat')
xlabel('x'); ylabel( 'phi');
title('final no norm')
colorbar
axis square

subplot(2,2,3)
imagesc( x, phi, ginhom')
xlabel('x'); ylabel( 'phi');
title('inhom')
colorbar
axis square

subplot(2,2,4)
imagesc( x, phi, ghom')
xlabel('x'); ylabel( 'phi');
title('hom')
colorbar
axis square

legend( 'original', 'perturbation', 'final')

%% 3D

% Guassian maker with center as input

