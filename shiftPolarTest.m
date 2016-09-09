% Add Subroutine path
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );

run1D = 0;
run2D = 0;
run3D = 1;

if run1D
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
pVar = 2 * pi / 4;

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
end

%% 2D example'

if run2D

% Parameters
Nx  = 24;
Nm  = 24;
Lx  = 1;
spatC = 2;

%Perturbation parameters
aXf = 0.5; % factor of ampitude you want pert to be
aPhif = 0.5; % factor of ampitude you want pert to be
varX = Lx / 4;
varPhi = 2 * pi / 4;
centerX = Lx / 2;

% grid stuff
dphi = 2 * pi / Nm;
dx = Lx / Nx;
phi = 0 : dphi : 2*pi - dphi;
x = 0 : dx : Lx - dx;

% 1D distro
isoC = 1 / (2 .* pi);
f = 1 + cos( 2 .* phi );
normTemp = trapz_periodic( phi, f );
f =  1 ./ normTemp .* f;
%Perturbation parameters
aX = aXf .* spatC;
aPhi   = aPhif .* isoC; 
centerPhi = phi( f == max(f) );
if length(centerPhi) > 1; centerPhi = centerPhi(1); end;

% Build 2D
f = repmat( f, [Nx,1] );
rho = spatC .* f;

% build perturbations
pX = gaussShift( x, aX, varX , centerX); % amp(x)
pX2 = repmat( pX', [1, Nm] );
pPhi = gaussShift( phi, aPhi, varPhi, centerPhi ); % amp(x)
pPhi2 = repmat( pPhi, [Nx, 1] );

pPolarAll    = pPhi2;
pPolarStripe = pX2 .* pPhi2;

% Homogenous polar everywhere
gHomPAll = rho + pPolarAll;
minG = min( min( gHomPAll ) );
if minG < 0;
  gHomPAll = gHomPAll - minG;
end
normTemp = trapz_periodic( x, trapz_periodic( phi, gHomPAll, 2 ) );
gHomPAll = spatC ./ normTemp .* gHomPAll;

% Homogenous polar stripe
gHomPStr = rho + pPolarStripe;
minG = min( min( gHomPStr ) );
if minG < 0;
  gHomPStr = gHomPStr - minG;
end
% Normalize the concentration
normVec = trapz_periodic( phi, gHomPStr, 2 );
normRep = repmat( normVec, [1, Nm ] );
gHomPStr = spatC ./ normRep .* gHomPStr;

% InHomogenous polar stripe
gInhomPStr = rho + pPolarStripe;
minG = min( min( gInhomPStr ) );
if minG < 0;
  gInhomPStr = gInhomPStr - minG;
end
% Normalize the concentration
normTemp = trapz_periodic( x, trapz_periodic( phi, gInhomPStr, 2 ), 1 );
gInhomPStr = spatC ./ normTemp .* gInhomPStr;
% % homogenous

figure()
subplot(2,3,1)
imagesc( x, phi, rho')
xlabel('x'); ylabel( 'phi');
title( 'orig' )
colorbar
axis square

subplot(2,3,2)
imagesc( x, phi, gHomPAll')
xlabel('x'); ylabel( 'phi');
title('Polar everywhere')
colorbar
axis square

subplot(2,3,3)
imagesc( x, phi, gHomPStr')
xlabel('x'); ylabel( 'phi');
title('Homogenous Stripe')
colorbar
axis square

subplot(2,3,4)
imagesc( x, phi, gInhomPStr')
xlabel('x'); ylabel( 'phi');
title('Inhomogenous Stripe')
colorbar
axis square


subplot(2,3,5)
imagesc( x, phi, pPolarAll')
xlabel('x'); ylabel( 'phi');
title('Polar All pert')
colorbar
axis square

subplot(2,3,6)
imagesc( x, phi, pPolarStripe')
xlabel('x'); ylabel( 'phi');
title('Inhom Stripe perturb')
colorbar
axis square

%%
figure()
plot( phi, gHomPAll( 1, : ), phi, gHomPAll( Nx/4, : ), phi, gHomPAll( Nx/2, : ) )
title('gHomPolarAll')

figure()
plot( phi, gHomPStr( 1, : ), phi, gHomPStr( Nx/4, : ), phi, gHomPStr( Nx/2, : ) )
title('gHom Stripe')
%%
end
%% 3D
if run3D
% Parameters
Nx  = 24;
Ny  = 24;
Nm  = 24;
Lx  = 1;
Ly  = 1;
spatC = 2;

%Perturbation parameters
aXf = 0.5; % factor of ampitude you want pert to be
aYf = 0.5;
aPhif = 0.5; % factor of ampitude you want pert to be
varX = Lx / 4;
varY = Ly / 4;
centerX = Lx / 2;
centerY = Ly / 2;
varPhi = 2 * pi / 4;

% grid stuff
dphi = 2 * pi / Nm;
phi = 0 : dphi : 2*pi - dphi;
dx = Lx / Nx;
x = 0 : dx : Lx - dx;
dy = Ly / Ny;
y = 0 : dy : Ly - dy;

% 1D distro
isoC = 1 / (2 .* pi);
f = 1 + cos( 2 .* phi );
normTemp = trapz_periodic( phi, f );
f =  1 ./ normTemp .* f;
%Perturbation parameters
aX = aXf .* spatC;
aY = aYf .* spatC;
aPhi   = aPhif .* isoC; 
centerPhi = phi( f == max(f) );
if length(centerPhi) > 1; centerPhi = centerPhi(1); end;

% Build 2D
f = repmat( reshape(f, [1,1,Nm] ), [Nx,Ny,1] );
rho = spatC .* f;

% build perturbations
pX = gaussShift( x, aX, varX , centerX); % amp(x)
pX3 = repmat( pX', [1, Ny, Nm] );
pY = gaussShift( y, aY, varY , centerY); % amp(x)
pY3 = repmat( pY, [Nx, 1, Nm] );
pPhi = gaussShift( phi, aPhi, varPhi, centerPhi ); % amp(x)
pPhi3 = repmat( reshape( pPhi, [1, 1, Nm] ), [Nx, Ny, 1] );

pPolarAll    = pPhi3;
pPolarStripe = pX3 .* pPhi3;
pPolarBall = pX3 .* pY3 .* pPhi3;

% Homogenous polar everywhere
gHomPAll = rho + pPolarAll;
minG = min( min( min( gHomPAll ) ) );
if minG < 0;
  gHomPAll = gHomPAll - minG;
end
normTemp = trapz_periodic( x , trapz_periodic( y, trapz_periodic( phi, gHomPAll, 3 ), 2 ), 1 );
gHomPAll = spatC ./ normTemp .* gHomPAll;

% Homogenous polar stripe
gHomPStr = rho + pPolarStripe;
minG = min( min( min( gHomPStr ) ) );
if minG < 0;
  gHomPStr = gHomPStr - minG;
end
% Normalize the concentration
normSqr = trapz_periodic( phi, gHomPStr, 3 );
normRep = repmat( normSqr, [1, 1, Nm ] );
gHomPStr = spatC ./ normRep .* gHomPStr;

% InHomogenous polar stripe
gInhomPStr = rho + pPolarStripe;
minG = min ( min( min( gInhomPStr ) ) );
if minG < 0;
  gInhomPStr = gInhomPStr - minG;
end
% Normalize the concentration
normTemp = trapz_periodic( x , trapz_periodic( ...
  y, trapz_periodic( phi, gInhomPStr, 3 ), 2 ), 1 );
gInhomPStr = spatC ./ normTemp .* gInhomPStr;

% Homogenous polar ball
gHomPBall = rho + pPolarBall;
minG = min( min( min( gHomPBall ) ) );
if minG < 0;
  gHomPBall = gHomPBall - minG;
end
% Normalize the concentration
normSqr = trapz_periodic( phi, gHomPBall, 3 );
normRep = repmat( normSqr, [1, 1, Nm ] );
gHomPBall = spatC ./ normRep .* gHomPBall;

% InHomogenous polar ball
gInhomPBall = rho + pPolarBall;
minG = min ( min( min( gInhomPBall ) ) );
if minG < 0;
  gInhomPBall = gInhomPBall - minG;
end
% Normalize the concentration
normTemp = trapz_periodic( x , trapz_periodic( ...
  y, trapz_periodic( phi, gInhomPBall, 3 ), 2 ), 1 );
gInhomPBall = spatC ./ normTemp .* gInhomPBall;


%% figures
gHomPAll; gHomPStr; gInhomPStr; gHomPBall; gInhomPBall; 
gTemp = gInhomPBall;
gStr = 'inhomo polar ball';

[~,~,phi3D] = meshgrid(phi);
cosPhi3d = cos(phi3D);
sinPhi3d = sin(phi3D);
cos2Phi3d = cosPhi3d .^ 2;
sin2Phi3d = sinPhi3d .^ 2;
cossinPhi3d = cosPhi3d .* sinPhi3d;

[C,PO,~,~,NO,~,~] = OpCPNCalc(Nx,Ny,gTemp,...
phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);


figure()
subplot(2,2,1)
imagesc( x, y, C')
xlabel('x'); ylabel( 'y');
title( 'C' )
colorbar
axis square

subplot(2,2,2)
imagesc( x, y, PO')
xlabel('x'); ylabel( 'y')
title('PO')
ax = gca;
ax.CLim = [0 max(max( PO ) ) ];
colorbar
axis square


subplot(2,2,3)
imagesc( x, y, NO')
xlabel('x'); ylabel( 'y');
ax = gca;
ax.CLim = [0 max(max( NO ) ) ];
title('NO')
colorbar
axis square

subplot(2,2,4)
plot( phi, reshape( gTemp( 1, 1, : ), [Nm, 1] ), ...
  phi, reshape( gTemp( 1, Ny/2, : ), [Nm, 1] ), ...
  phi, reshape( gTemp( Nx/2, 1, : ), [Nm, 1] ), ...
  phi, reshape( gTemp( Nx/2, Ny/2, : ), [Nm, 1] ) )
xlabel('phi'); ylabel( 'distro');
title( gStr );
legend( ' (1,1) ', ' (1,Ny/2) ', ' Nx/2,1 ', ' Nx/2, Ny/2 ' );
axis square

 
end
