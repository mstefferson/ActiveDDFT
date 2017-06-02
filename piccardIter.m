% System parameters
n = 64;
l = 10;
c = 2;
alpha = 0.5; %mixing parameter
maxIters = 50; %mixing parameter
randEps = 0.001;
ssEps = 0.00001;
lrLs1 = 1;
lrLs2 = 1.855;
lrEs1 = 1;
lrEs2 = 1;
% set other parameters
n1 = n;
n2 = n;
l1 = l;
l2 = l;
convolFact = l1 * l2 / (n1 * n2 );
% paths
addpath( genpath( './src') );
% initialize rho
rho = c * ones( n1, n2);
rhoBulk = rho;
rho = rho + randEps * rand(n1,n2);
rhoInit = rho;
dRho = 100;
% get potential
[~, vFt] = softShoulder2D( lrLs1, lrLs2, lrLs1, lrLs2, n1, l1, n2, l2 );
% interate to steady state
counter = 0;
while dRho > ssEps
  % calculate correlation
  cFt = -convolFact .* vFt .* fftshift( fftn( rho ) );
  c = real( ifftn( ifftshift( cFt ) ) );
  % implement piccard iteration
  rhoTilde = rhoBulk .* exp( c );
  rhoNext = (1-alpha) .* rhoTilde + alpha * rho;
  % update
  dRho = max( abs( rho(:) - rhoNext(:) ) );
  rho = rhoNext;
  counter = counter + 1;
  if counter >= maxIters
    fprintf('Reached max iterations %d\n', maxIters);
  end
end
% display info
fprintf('Ran %d interactions. dRho = %f \n', counter, dRho);
figure()
subplot(1,2,1)
imagesc( rhoInit );
colorbar

subplot(1,2,2)
imagesc( rho);
colorbar

