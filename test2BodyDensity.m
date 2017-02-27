currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% system obj
systemObj.n1 = 16;
systemObj.n2 = 16;
systemObj.n3 = 16;
systemObj.l1 = 3;
systemObj.l2 = 3;
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
pCorr1 = 1 + 1 / (2 * pi) * mayerInt;
% now with a density
%load( [path '/op_' filename '.mat'] )
%% version 2
% rhoTemp = rand( systemObj.n1, systemObj.n2, systemObj.n3 );
rhoTemp = 3 .* ones( systemObj.n1, systemObj.n2, systemObj.n3 );
%load('rhoTest30.mat');
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
pCorr2 = ones( systemObj.n1, systemObj.n2,systemObj.n1, systemObj.n2);
%% version 1 wrong phi integral
for ii = 1:systemObj.n1
  % inds for loop over delta x
  inds1 = mod( (ii-1) + indsDelta1 - 1, systemObj.n1 ) + 1;
  for jj = 1:systemObj.n2
    % inds for loop over delta y
    inds2 = mod( (jj-1) + indsDelta2 - 1, systemObj.n2 ) + 1;
    for kk = 1:totalInds1
      ii2 = inds1(kk);
      dii = indsDelta1(kk);
      for ll = 1:totalInds2
        % plug into
        jj2 = inds2(ll);
        djj = indsDelta2(ll);
        vec2Int = rhoTemp(ii, jj, :) .* rhoTemp( ii2, jj2, :)...
          .* mayer( dii, djj, : );
        vec2Int = reshape( vec2Int, [1, systemObj.n3] );
        num = trapz_periodic( phi, vec2Int );
        den = ( cTemp(ii,jj) .* cTemp(ii2, jj2) );
        pCorr2( ii, jj, ii2, jj2 ) = pCorr2( ii, jj, ii2, jj2 ) + num ./ den;
      end
    end
  end
end

%% version 2 wrong phi integral
pCorr2v2 = ones( systemObj.n1, systemObj.n2,systemObj.n1, systemObj.n2);
for ii = 1:systemObj.n1
  % inds for loop over delta x
  inds1 = mod( (ii-1) + indsDelta1 - 1, systemObj.n1 ) + 1;
  for jj = 1:systemObj.n2
    % inds for loop over delta y
    inds2 = mod( (jj-1) + indsDelta2 - 1, systemObj.n2 ) + 1;
    vec2Int = rhoTemp(ii, jj, :) .* rhoTemp( inds1, inds2, :)...
      .* mayer( indsDelta1, indsDelta2, : );
    int2 = trapz_periodic( phi, mayer( indsDelta1, indsDelta2, : ), 3 );
%     vec2Int = reshape( vec2Int, [totalInds1, totalInds2, systemObj.n3] );
    num = trapz_periodic( phi, vec2Int, 3 );
    den = ( cTemp(ii,jj) .* cTemp(inds1, inds2) );
    newterm = num ./ den;
    reshapeNewterm = reshape( newterm, [1, 1, 13, 13] );
    pCorr2v2( ii, jj, inds1, inds2 ) = pCorr2v2( ii, jj, inds1, inds2 ) + reshapeNewterm;
  end
end

pCorr2correcInt = ones( systemObj.n1, systemObj.n2,systemObj.n1, systemObj.n2);
%% version 1 wrong phi integral
indsDelta3 = 1:systemObj.n3;
vec2Int = zeros( 1, systemObj.n3 );
for ii = 1:systemObj.n1
  % inds for loop over delta x
  inds1 = mod( (ii-1) + indsDelta1 - 1, systemObj.n1 ) + 1;
  for jj = 1:systemObj.n2
    % inds for loop over delta y
    inds2 = mod( (jj-1) + indsDelta2 - 1, systemObj.n2 ) + 1;
    for kk = 1:totalInds1
      ii2 = inds1(kk);
      dii = indsDelta1(kk);
      for ll = 1:totalInds2
        % plug into
        jj2 = inds2(ll);
        djj = indsDelta2(ll);
        for mm = 1:systemObj.n3
          inds3 = mod( (mm-1) + indsDelta3 - 1, systemObj.n3 ) + 1;
          vec2Temp = rhoTemp(ii, jj, mm) .* rhoTemp( ii2, jj2, inds3)...
            .* mayer( dii, djj, indsDelta3 );
          vec2Temp = reshape( vec2Temp, [1, systemObj.n3] );
          vec2Int(mm) = trapz_periodic( phi, vec2Temp );
        end
        num = trapz_periodic( phi, vec2Int );
        den = ( cTemp(ii,jj) .* cTemp(ii2, jj2) );
        pCorr2correcInt( ii, jj, ii2, jj2 ) = pCorr2correcInt( ii, jj, ii2, jj2 ) + num ./ den;
      end
    end
  end
end
%% version 2 correct phi integral
pCorr2v2correctInt = ones( systemObj.n1, systemObj.n2,systemObj.n1, systemObj.n2);
vec2Int = zeros( totalInds1, totalInds2, systemObj.n3 );
for ii = 1:systemObj.n1
  % inds for loop over delta x
  inds1 = mod( (ii-1) + indsDelta1 - 1, systemObj.n1 ) + 1;
  for jj = 1:systemObj.n2
    % inds for loop over delta y
    inds2 = mod( (jj-1) + indsDelta2 - 1, systemObj.n2 ) + 1;
    for mm = 1:systemObj.n3
      inds3 = mod( (mm-1) + indsDelta3 - 1, systemObj.n3 ) + 1;
      vec2Temp = rhoTemp(ii, jj, mm) .* rhoTemp( inds1, inds2, inds3)...
        .* mayer( indsDelta1, indsDelta2, indsDelta3 );
      vec2Int(:,:,mm) = trapz_periodic( phi, vec2Temp, 3 );
    end
%     vec2Int = reshape( vec2Int, [totalInds1, totalInds2, systemObj.n3] );
    num = trapz_periodic( phi, vec2Int, 3 );
    den = ( cTemp(ii,jj) .* cTemp(inds1, inds2) );
    newterm = num ./ den;
    reshapeNewterm = reshape( newterm, [1, 1, 13, 13] );
    pCorr2v2correctInt( ii, jj, inds1, inds2 ) = pCorr2v2correctInt( ii, jj, inds1, inds2 ) + reshapeNewterm;
  end
end
%% Average
aveCorr2v1 = zeros( systemObj.n1, systemObj.n2);
n1n2 = systemObj.n1 .* systemObj.n2;
for ii = indsDelta1
  %ind2Ave1 = mod( allInds1 + (ii-1), systemObj.n1 );
  for jj = 1:indsDelta2
    %ind2Ave2 = mod( allInds2 + (jj-1), systemObj.n2 );
    for kk = 1:systemObj.n1
      ii2 = mod( kk + (ii-2), systemObj.n1 ) + 1;
      for ll = 1:systemObj.n2
        jj2 = mod( ll + (jj-1), systemObj.n2 ) + 1;
        %aveCorrelation2(ii,jj) = aveCorrelation2(ii,jj) + ...
        %fullCorrelation2( ind2Ave1(kk), inds2Ave2(ll), allInds1(ii), allInds2(jj)  ) ;
        aveCorr2v1(ii,jj) = aveCorr2v1(ii,jj) + ...
          pCorr2( kk, ll, ii2, jj2  ) ;
      end
    end
    aveCorr2v1(ii,jj) = aveCorr2v1(ii,jj) ./ ( n1n2 );
  end
end
%% Average 1.5
aveCorr2v12 = zeros( systemObj.n1, systemObj.n2);
kk = 1:systemObj.n1;
ll = 1:systemObj.n2;
for ii = indsDelta1
  %ind2Ave1 = mod( allInds1 + (ii-1), systemObj.n1 );
  for jj = 1:indsDelta2
    %ind2Ave2 = mod( allInds2 + (jj-1), systemObj.n2 );
    ii2 = mod( kk + (ii-2), systemObj.n1 ) + 1;
    jj2 = mod( ll + (jj-2), systemObj.n2 ) + 1;
    aveCorr2v12(ii,jj) = mean(mean(mean(mean( pCorr2( kk, ll, ii2, jj2  ) ))));
  end
end
%% Average version 2
allInds1 = 1:systemObj.n1;
allInds2 = 1:systemObj.n2;
aveCorr2v2 = zeros( systemObj.n1, systemObj.n2);
n1n2 = systemObj.n1 .* systemObj.n2;
for ii = 1:totalInds1
  ind2Ave1 = mod( allInds1 + (ii-2), systemObj.n1 ) + 1 ;
  for jj = 1:totalInds2
    ind2Ave2 = mod( allInds2 + (jj-2), systemObj.n2 ) + 1 ;
    %    matTemp = reshape( fullCorrelation2( allInds2, allInds2, ind2Ave1, ind2Ave2 ), [systemObj.n1, systemObj.n2] );
    tempMean = mean( mean( mean( mean( pCorr2( allInds2, allInds2, ind2Ave1, ind2Ave2 ) ) ) ) );
    %aveCorrelation2_v2(ii,jj) =  mean( mean( matTemp ) );
    aveCorr2v2(ii,jj) =  tempMean;
  end
end
%% center
center1 = round( systemObj.n2 / 2 ) + 1 ;
center2 = round( systemObj.n2 / 2 ) + 1 ;
pCorr1center = circshift( circshift( pCorr1, center1-1, 1 ), ...
  center2-1, 2 );
aveCorr2v1center = circshift( circshift( aveCorr2v1, center1-1, 1), ...
  center2-1, 2);
aveCorr2v12center = circshift( circshift( aveCorr2v12, center1-1, 1), ...
  center2-1, 2);
aveCorr2v2center = circshift( circshift( aveCorr2v2, center1-1, 1), ...
  center2-1, 2);
%% plot
%
figure(1)
imagesc( pCorr1center )
title('no rho dependence')
colorbar
%
figure(2)
imagesc( pCorr2(:,:,center1,center2) )
title('pair Correlations at fixed x2, y2')
colorbar
%
figure(3)
imagesc( pCorr2v2(:,:,center1,center2) )
title('pair Correlations at fixed x2, y2 version 2')
colorbar
%
figure(4)
imagesc( pCorr2correcInt(:,:,center1,center2) )
title('pair Correlations at fixed x2, y2 correct int')
colorbar
%
figure(5)
imagesc( pCorr2v2correctInt(:,:,center1,center2) )
title('pair Correlations at fixed x2, y2 version 2 correct int')
colorbar
% %
% figure(4)
% imagesc( aveCorr2v1center )
% title('rho dependence ave 1')
% colorbar
% %
% figure(5)
% imagesc( aveCorr2v12center )
% title('rho dependence ave 1.5')
% colorbar
% %
% figure(6)
% imagesc( aveCorr2v2center )
% title('rho dependence ave 2')
% colorbar
