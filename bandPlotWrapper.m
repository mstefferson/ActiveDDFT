%
addpath( genpath( './src' ) );
saveMe = 1;
plotMyBand = 1;
plotMySlice = 1;
plotMyPhase = 1;
dir2check = 'analyzedfiles/BandAnalysis';
if exist( 'bandTable', 'var') == 0
  [bandSum, bandTable] = bandAnalysis(dir2check);
end
% plot band
if plotMyBand
  % svg and pdfs for the band are crashing matlab...
  saveExtBand = 'png';
  filename = ...
    'Hr_rods_mayer_diag1_N160160160_ls1010_bc1.55_fD10_ICload_SM6_t20180222.02';
  fprintf('Starting plotBand\n')
  plotBand( [dir2check '/' filename] );
  fprintf('Finished plotBand\n')
  if saveMe
    saveTitle = 'band_example';   
    saveas( gcf, saveTitle, saveExtBand  )
    fprintf('Saved\n')
  end
end
% plot band slice
if plotMySlice
  saveExt = 'png';
  cWant = 1.55;
  fWant = [0, 5, 10 30];
  fprintf('Starting plotBandSlice\n')
  plotBandSlice( cWant, fWant, bandTable );
  fprintf('Finished plotBandSlice\n')
  if saveMe
    saveTitle = 'band_analysis';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end
% plot phase diagram
if plotMyPhase
  saveExt = 'svg';
  plotScaled = 0;
  plotTheory = 0;
  plotIN = 0;
  fprintf('Starting plotBandPhase\n')
  plotBandPhase( bandTable, plotIN, plotTheory, plotScaled )
  fprintf('Finished plotBandPhase\n')
  if saveMe
    saveTitle = 'phase_diagram';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end

