%
saveMe = 1;
plotMySlice = 1;
plotMyPhase = 1;
plotMyBand = 1;
saveExt = 'png';

dir2check = 'analyzedfiles/BandAnalysis';
if exist( 'bandTable', 'var') == 0
  [bandSum, bandTable] = bandAnalysis(dir2check);
end
% plot band slice
if plotMySlice
  cWant = 1.55;
  fWant = [0, 5, 10 30];
  plotBandSlice( cWant, fWant, bandTable );
  if saveMe
    saveTitle = 'band_analysis';
    saveas( gcf, saveTitle, saveExt )
  end
end
% plot phase diagram
if plotMyPhase
  plotScaled = 0;
  plotTheory = 0;
  plotBandPhase( bandTable, plotTheory, plotScaled )
  if saveMe
    saveTitle = 'phase_diagram';
    saveas( gcf, saveTitle, saveExt )
  end
end
if plotMyBand
  % plot band
  filename = ...
    'Hr_rods_mayer_diag1_N160160160_ls1010_bc1.55_fD10_ICload_SM6_t20180222.02';
  plotBand( [dir2check '/' filename] );
  if saveMe
    saveTitle = 'band_example';
    saveas( gcf, saveTitle, saveExt )
  end
end

