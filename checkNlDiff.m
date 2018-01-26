function checkNlDiff( n3, gridObj, GammaCube_FT, lop, rho_FT, rho, ...
  polarDrive, densityDepDiff )
  gammaExReal = real( ifftn( ifftshift( GammaCube_FT ) ) );
  gammaExSqReal = trapz_periodic( gammaExReal, 3 );
  gammaExPhiReal =  reshape( trapz_periodic( trapz_periodic(...
    gammaExReal,1 ),2 ), [1 n3] );
  gammaDiff = lop .* rho_FT;
  gammaDiffReal = real( ifftn( ifftshift( gammaDiff ) ) );
  gammaDiffSqReal =  trapz_periodic( gammaDiffReal, 3 );
  gammaDiffPhiReal =  reshape( trapz_periodic( trapz_periodic(...
    gammaDiffReal,1 ),2 ), [1 n3] );
  c = trapz_periodic( rho, 3 );
  f = reshape( trapz_periodic( trapz_periodic( rho, 1 ), 2 ), [1 n3] );
  % plot it
  %% spatial
  figure()
  subplot(2,2,1)
  imagesc( c )
  colorbar
  title('c')
  subplot(2,2,2)
  imagesc( gammaDiffSqReal )
  title('Diffusion')
  colorbar
  subplot(2,2,3)
  imagesc( gammaExSqReal )
  title('Excess')
  colorbar
  subplot(2,2,4)
  sumGamm = gammaDiffSqReal + gammaExSqReal;
  imagesc( sumGamm )
  title('Sum')
  colorbar
  %% angle
  figure()
  subplot(2,1,1)
  plot( gridObj.x3, gammaDiffPhiReal,...
    gridObj.x3, gammaExPhiReal,...
    gridObj.x3, gammaExPhiReal + gammaDiffPhiReal)
  ax = gca;
  ax.XLim = [gridObj.x3(1) gridObj.x3(end)];
  legend('Diffusion', 'Excess', 'Sum')
  title( ['n = ' num2str(n3,'%d')] )
  subplot(2,1,2)
  plot(gridObj.x3,f)
  ax = gca;
  ax.XLim = [gridObj.x3(1) gridObj.x3(end)];
  title('f')
  %% plot slices
  figure
  n1Center = gridObj.k1ind0;
  n2Center = gridObj.k2ind0; 
  subplot(2,2,1)
  plot( gridObj.x1, c(:,n2Center), gridObj.x2, c(n1Center,:) )
  legend('along n1', 'along n2')
  title('c')
  subplot(2,2,2)
  plot( gridObj.x1, gammaDiffSqReal(:,n2Center), gridObj.x2, gammaDiffSqReal(n1Center,:) )
  legend('along n1', 'along n2')
  title('Diffusion')
  subplot(2,2,3)
  plot( gridObj.x1, gammaExSqReal(:,n2Center), gridObj.x2, gammaExSqReal(n1Center,:) )
  legend('along n1', 'along n2')
  title('Excess')
  subplot(2,2,4)
  sumGamm = gammaDiffSqReal + gammaExSqReal;
  plot( gridObj.x1, sumGamm(:,n2Center), gridObj.x2, sumGamm(n1Center,:) )
  legend('along n1', 'along n2')
  title('Sum')
  %% Fluxes (all angles should be the same if no angular dep in rho
  phiInd = 50;
  iotaDiff1 = densityDepDiff.calcIotaDiff( rho_FT, densityDepDiff.Ik{1} );
  iotaDiff2 = densityDepDiff.calcIotaDiff( rho_FT, densityDepDiff.Ik{2} );
  figure()
  subplot(2,3,1)
  imagesc( polarDrive.Iota1(:,:,phiInd) );
  colorbar
  title('Drive 1')
  subplot(2,3,2)
  imagesc( iotaDiff1(:,:,phiInd));
  colorbar
  title('Diff 1')
  subplot(2,3,3)
  sumIota1 = iotaDiff1(:,:,phiInd) + polarDrive.Iota1(:,:,phiInd);
  imagesc( sumIota1 );
  colorbar
  title('Total 1')
  subplot(2,3,4)
  imagesc( polarDrive.Iota2(:,:,phiInd) )
  title('Drive 2')
  colorbar
  subplot(2,3,5)
  imagesc( iotaDiff2(:,:,phiInd) )
  title('Diff 2')
  colorbar
  subplot(2,3,6)
  sumIota2 = iotaDiff2(:,:,phiInd) + polarDrive.Iota2(:,:,phiInd);
  imagesc( sumIota2 );
  colorbar
  title('Total 2')
  
  %%
  figure()
  imagesc( densityDepDiff.DNlPos ); colorbar
%     subplot(1,3,3)
%   imagesc( densityDepDiff.DNlPos )
  %%
  keyboard
