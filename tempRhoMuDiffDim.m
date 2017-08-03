%%
fprintf('\n\n Angular 1D potential\n\n')
n1 = 100;
n2 = 100;
n3 = 100;
nt = 10;
rho = rand(n1,n2,n3);
rhoFt = fftshift( fftn( rho ) );
[n3mesh] = 1:n3;
v1 = exp( -n3mesh .^2 );

%% way 1
fprintf('Straight Multiply mu\n')
ind1 = n1/2 + 1;
ind2 = n2/2 + 1;
ind3 = 1:n3;
reshapeInds = [1 1 n3];
vFt = fftshift(fftn(v1));
vFt = reshape( vFt, reshapeInds );
tic
for ii = 1:nt
  muFt = 1 / (n1*n2*n3) * rhoFt(ind1,ind2,ind3) .* vFt;
  mu = real( ifftn( ifftshift( muFt ) ) );
  j1 = rho .* mu;
end
toc
%%

%% way 2
fprintf('Rep mu to missing dimension\n')
tic
reshapeInds = [1 1 n3];
repInds = [n1 n2 1];
ind1 = n1/2 + 1;
ind2 = n2/2 + 1;
ind3 = 1:n3;
vFt = fftshift(fftn(v1));
vFt = reshape( vFt, reshapeInds );
for ii = 1:nt
  muFt = 1 / (n1*n2*n3) * rhoFt(ind1,ind2,ind3) .* vFt;
  mu = real( ifftn( ifftshift( muFt ) ) );
  mu = repmat( mu, repInds );
  j2 = rho .* mu;
end
toc

%% way 3
fprintf('Build and fill vFt to missing dimension\n')
tic
vFt = fftshift(fftn(v1));
ind1 = n1/2 + 1;
ind2 = n2/2 + 1;
ind3 = 1:n3;
vFtRep = zeros( n1,n2,n3 );
vFtRep( ind1,ind2,ind3) = (n1 *n2) * vFt;
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFtRep;
  mu = real( ifftn( ifftshift( muFt ) ) );
  j3 = rho .* mu;
end
toc

%% way 4
fprintf('Rep v to missing dimension\n')
tic
repInds = [n1,n2,1];
v1Reshape = reshape( v1, reshapeInds );
vRep = repmat(v1Reshape, repInds);
vFtRep = fftshift(fftn(vRep));
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFtRep;
  mu = real( ifftn( ifftshift( muFt ) ) );
  j4 = rho .* mu;
end
toc
%%
fprintf('checking: %f, %f, %f \n', ...
  max( abs( j1(:) - j2(:) ) ), ...
  max( abs( j1(:) - j3(:) ) ), ...
  max( abs( j1(:) - j4(:) ) ));

%%%%%%% 2D Pot %%%%%%%
%%
fprintf('\n\n Spatial 2D potential\n\n')
rho = rand(n1,n2,n3);
rhoFt = fftshift( fftn( rho ) );
[n1mesh, n2mesh] = meshgrid( 1/(n3/10)*(1:n3) );
v2 = exp( -n1mesh .^2 -n2mesh .^2 );

%% way 1
fprintf('Straight Multiply mu\n')
tic
ind1 = 1:n1;
ind2 = 1:n1;
ind3 = n3/2+1;
reshapeInds = [n1, n2, 1];
vFt = fftshift(fftn(v2));
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt(ind1,ind2,ind3) .* vFt;
  mu1 = real(  ifftn( ifftshift( muFt ) ) );
  mu1 = reshape( mu1, reshapeInds );
  j1 = rho .* mu1;
end
toc
%%

%% way 2
fprintf('Rep mu to missing dimension\n')
tic
ind1 = 1:n1;
ind2 = 1:n1;
ind3 = n3/2+1;
repInds = [1,1,n3];
vFt = fftshift(fftn(v2));
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt(:,:,n3/2+1) .* vFt;
  mu2 = real( ifftn( ifftshift( muFt ) ) );
  mu2 = repmat( mu2, repInds );
  j2 = rho .* mu2;
end
toc

%% way 3
fprintf('Build and fill vFt to missing dimension\n')
tic
vFt = fftshift(fftn(v2));
vFtRep = zeros( n1,n2,n3 );
vFtRep( :,:,n3/2+1 ) = n3 * vFt;
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFtRep;
  mu3 = real(  ifftn( ifftshift( muFt ) ) );
  j3 = rho .* mu3;
end
toc

%% way 4
fprintf('Rep v to missing dimension\n')
tic
repInds = [1,1,n3];
vRep = repmat(v2, repInds);
vFtRep = fftshift(fftn(vRep));
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFtRep;
  mu4 = real( ifftn( ifftshift( muFt ) ) );
  j4 = rho .* mu4;
end
toc

fprintf('checking: %f, %f, %f \n', ...
  max( abs( j1(:) - j2(:) ) ), ...
  max( abs( j1(:) - j3(:) ) ), ...
  max( abs( j1(:) - j4(:) ) ));

%%%%%%% Multiple potentials %%%%%%%%%

fprintf('\n\n Mult potentials \n\n')
rho = rand(n1,n2,n3);
rhoFt = fftshift( fftn( rho ) );
[n1mesh, n2mesh] = meshgrid( 1/(n3/10)*(1:n3) );
v2 = exp( -n1mesh .^2 -n2mesh .^2 );
[n1mesh, n2mesh, n3mesh] = meshgrid( 1/(n3/10)*(1:n3) );
v3 = exp( -n1mesh .^2 -n2mesh .^2 -n3mesh .^2 );

%% way 1
fprintf('Straight Multiply mu\n')
tic
v2Ft = fftshift(fftn(v2));
v3Ft = fftshift(fftn(v3));
for ii = 1:nt
  mu2Ft = 1 / (n1*n3*n3) * rhoFt(:,:,n3/2+1) .* v2Ft;
  mu2 = real( ifftn( ifftshift( mu2Ft ) ) );
  mu2 = reshape( mu2, [n1,n2,1] );
  mu3Ft = 1 / (n1*n3*n3) * rhoFt .* v3Ft;
  mu3 = real( ifftn( ifftshift( mu3Ft ) ) );
  mu = mu2 + mu3;
  j1 = rho .* mu;
end
toc
%%

%% way 2
fprintf('Rep mu to missing dimension\n')
tic
v2Ft = fftshift(fftn(v2));
repInds = [1,1,n3];
for ii = 1:nt
  mu2Ft = 1 / (n1*n3*n3) * rhoFt(:,:,n3/2+1) .* v2Ft;
  mu2 = real( ifftn( ifftshift( mu2Ft ) ) );
  mu2 = repmat( mu2, repInds );
  mu3Ft = 1 / (n1*n3*n3) * rhoFt .* v3Ft;
  mu3 = real( ifftn( ifftshift( mu3Ft ) ) );
  mu = mu2 + mu3;
  j2 = rho .* mu;
end
toc

%% way 3
fprintf('Build and fill vft2 and combine vFt \n')
tic
v2Ft = fftshift(fftn(v2));
v2FtRep = zeros( n1,n2,n3 );
v2FtRep( :,:,n3/2+1 ) = n3 * v2Ft;
vFt = v2FtRep + v3Ft;
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFt;
  mu = real( ifftn( ifftshift( muFt ) ) );
  j3 = rho .* mu;
end
toc

%% way 4
fprintf('Combine v straight add\n')
tic
v2reshape = reshape(v2, [n1, n2, 1] );
v = v2reshape + v3;
vFt = fftshift(fftn(v));
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFt;
  mu = real( ifftn( ifftshift( muFt ) ) );
  j4 = rho .* mu;
end
toc

%% way 5
fprintf('Combine v rep add\n')
tic
v2rep = repmat( v2, [1, 1, n3] );
v = v2rep + v3;
vFt = fftshift(fftn(v));
for ii = 1:nt
  muFt = 1 / (n1*n3*n3) * rhoFt .* vFt;
  mu = real( ifftn( ifftshift( muFt ) ) );
  j5 = rho .* mu;
end
toc

fprintf('checking: %f, %f, %f, %f \n', ...
  max( abs( j1(:) - j2(:) ) ), ...
  max( abs( j1(:) - j3(:) ) ), ...
  max( abs( j1(:) - j4(:) ) ), ...
  max( abs( j1(:) - j5(:) ) ));
