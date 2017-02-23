currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% system obj
systemObj.n1 = 32;
systemObj.n2 = 32;
systemObj.n3 = 32;
systemObj.l1 = 10;
systemObj.l2 = 10;
particleObj.lMaj = 1;
% grab files
files = {
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD0_IC1_SM6_t09.04';
  };
%filename = files{ii};
%path = ['./analyzedfiles/summarized/vDsweep/' filename ];
%load( [path '/params_' filename '.mat'] )
% build phi
phi = 0 : 2*pi / systemObj.n3 : 2*pi * ( 1 - 1 ./ systemObj.n3 );
% get MayerFunction
[mayer] = mayerFncHr(...
  systemObj.n1, systemObj.n2, systemObj.n3, ...
  systemObj.l1, systemObj.l2, particleObj.lMaj) ;
%% intergrate mayer
mayerInt = trapz_periodic( phi, mayer, 3 );
% non-density dependent correlation integrate g(r1,u1,r2,u2)
pairCorrelation1 = 4 * pi ^ 2 + 2 * pi * mayerInt;
% now with a density
%load( [path '/op_' filename '.mat'] )
%% version 2
rhoTemp = rand( systemObj.n1, systemObj.n2, systemObj.n3 );
load('rhoTest30.mat');
% rhoTemp = rho;
cTemp = trapz_periodic( phi, rhoTemp, 3 );
% cTemp = C_rec(:,:,end) ./ pi;
%rhoTemp = denRecObj.rhoFinal(:,:,:);
% try and integrate
maxDelta1 =  ceil( systemObj.n1 .* particleObj.lMaj ./ systemObj.l1 );
indsDelta1 = mod( -maxDelta1:(maxDelta1+1) - 1, systemObj.n1 ) + 1;
totalInds1 = 2 .* maxDelta1 + 1;
maxDelta2 =  ceil( systemObj.n2 .* particleObj.lMaj ./ systemObj.l2 );
indsDelta2 = mod( -maxDelta2:(maxDelta2+1) - 1, systemObj.n1 ) + 1;
totalInds2 = 2 .* maxDelta2 + 1;
fullCorrelation2 = zeros( systemObj.n1, systemObj.n2,systemObj.n1, systemObj.n2);
for ii = 1:systemObj.n1
  % inds for loop over delta x
  inds1 = mod( (ii-1) + indsDelta1 - 1, systemObj.n1 ) + 1;
  for jj = 1:systemObj.n2
    % inds for loop over delta y
    inds2 = mod( (jj-1) + indsDelta2 - 1, systemObj.n2 ) + 1;
    for kk = 1:totalInds1
      for ll = 1:totalInds2
        % plug into
        vec2Int = rhoTemp(ii, jj, :) .* rhoTemp( inds1(kk), inds2(ll), :)...
          .* mayer( indsDelta1(kk), indsDelta2(ll), : );
        vec2Int = reshape( vec2Int, [1, systemObj.n3] );
        fullCorrelation2( ii, jj, inds1(kk), inds2(ll) ) = 2 * pi * trapz_periodic( phi, vec2Int );
        fullCorrelation2( ii, jj, inds1(kk), inds2(ll) ) = fullCorrelation2( ii, jj, inds1(kk), inds2(ll) ) ./ ( cTemp(ii,jj) .* cTemp(inds1(kk), inds2(ll)) );
      end
    end
  end
end
%% Average
allInds1 = 1:systemObj.n1;
allInds2 = 1:systemObj.n2;
aveCorrelation2 = zeros( systemObj.n1, systemObj.n2);
n1n2 = systemObj.n1 .* systemObj.n2;
for ii = 1:totalInds1
  ind2Ave1 = mod( allInds1 + (ii-1), systemObj.n1 );
  for jj = 1:totalInds2
    ind2Ave2 = mod( allInds2 + (jj-1), systemObj.n2 );
    for kk = 1:systemObj.n1
      deltaK = mod( kk + (ii-1), systemObj.n1 ) + 1;
      for ll = 1:systemObj.n2
        deltaL = mod( ll + (jj-1), systemObj.n2 ) + 1;
        %aveCorrelation2(ii,jj) = aveCorrelation2(ii,jj) + ...
          %fullCorrelation2( ind2Ave1(kk), inds2Ave2(ll), allInds1(ii), allInds2(jj)  ) ;
        aveCorrelation2(ii,jj) = aveCorrelation2(ii,jj) + ...
          fullCorrelation2( kk, ll, deltaK, deltaL  ) ;
      end
    end
    aveCorrelation2(ii,jj) = aveCorrelation2(ii,jj) ./ ( n1n2 );
  end
end
aveCorrelation2 = aveCorrelation2 + 1;
%% Average version 2
allInds1 = 1:systemObj.n1;
allInds2 = 1:systemObj.n2;
aveCorrelation2_v2 = zeros( systemObj.n1, systemObj.n2);
n1n2 = systemObj.n1 .* systemObj.n2;
for ii = 1:totalInds1
  ind2Ave1 = mod( allInds1 + (ii-2), systemObj.n1 ) + 1 ;
  for jj = 1:totalInds2
    ind2Ave2 = mod( allInds2 + (jj-2), systemObj.n2 ) + 1 ;
    matTemp = reshape( fullCorrelation2( allInds2, allInds2, ind2Ave1, ind2Ave2 ), [systemObj.n1, systemObj.n2] );
    tempMean = mean( mean( mean( mean( fullCorrelation2( allInds2, allInds2, ind2Ave1, ind2Ave2 ) ) ) ) );
    aveCorrelation2_v2(ii,jj) =  mean( mean( matTemp ) );
  end
end

aveCorrelation2 = aveCorrelation2 + 1;

