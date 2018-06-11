%%
plotMyBand = 0;
plotMySlice = 0;
plotMySliceDeriv = 1;
plotMyPhase = 0;
plotMyIN = 0;
bandThres = 0.05;
saveMe = 0;
dockStyle = 'normal';
addpath( genpath( './src' ) );
dir2check = 'analyzedfiles/BandAnalysis';
overrideBandTable = 0;
if overrideBandTable
  fprintf('Overridding bandTable variable. Savin and using\n')
  [~, bandTable] = bandAnalysis(dir2check, bandThres);
  save('bandTable', 'bandTable')  
elseif exist( 'bandTable', 'var')
  fprintf('Found bandTable variable. Using that\n')
elseif exist( 'bandTable.mat', 'file')
  fprintf('Found saved bandTable file. Using that\n')
  load('bandTable.mat')
else
  fprintf('Cannot find bandTable variable. Saving and using\n')
  [~, bandTable] = bandAnalysis(dir2check, bandThres);
  save('bandTable', 'bandTable')
end

%%
% plot band
if plotMyBand
  % svg and pdfs for the band are crashing matlab...
  saveExtBand = 'png';
  filename = ...
    'Hr_rods_mayer_diag1_N160160160_ls1010_bc1.55_fD10_ICload_SM6_t20180222.02';
  fprintf('Starting plotBand\n')
  plotBand( [dir2check '/' filename] );
  fig = gcf;
  fig.WindowStyle = dockStyle;
  fprintf('Finished plotBand\n')
  if saveMe
    saveTitle = 'band_example';   
    saveas( gcf, saveTitle, saveExtBand  )
    fprintf('Saved\n')
  end
end
%% compare derivatives
if plotMySliceDeriv
  saveExt = 'png';
  cWant = 1.55;
  fWant = [0, 5, 10 30];
  fprintf('Starting plotBandSlice\n')
  plotCompareBandDeriv( cWant, fWant, bandTable );
  fprintf('Finished plotBandSlice\n')
  if saveMe
    saveTitle = 'band_analysis';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end
%% plot band slice
if plotMySlice
  saveExt = 'png';
  cWant = 1.55;
  fWant = [0, 5, 10 30];
  fprintf('Starting plotBandSlice\n')
  plotBandSlice( cWant, fWant, bandTable );
  fig = gcf;
  fig.WindowStyle = dockStyle;
  fprintf('Finished plotBandSlice\n')
  if saveMe
    saveTitle = 'band_analysis';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end

%% plot phase diagram
if plotMyPhase
  saveExt = 'svg';
  plotScaled = 0;
  plotTheory = 0;
  plotIN = 0;
  fprintf('Starting plotBandPhase\n')
  plotBandPhase( bandTable, 'cMax', plotIN, plotTheory )
  fig = gcf;
  fig.WindowStyle = dockStyle;
  fprintf('Finished plotBandPhase cPeak\n')
  if saveMe
    saveTitle = 'phase_diagram_peak';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
  plotBandPhase( bandTable, 'cFWHM', plotIN, plotTheory )
  fig = gcf;
  fig.WindowStyle = dockStyle;
  fprintf('Finished plotBandPhase cFWHM\n')
  if saveMe
    saveTitle = 'phase_diagram_width';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end

%% Phase sub
if plotMyPhase
  plotBandPhaseSub( bandTable, 'cMax','cFWHM', plotIN, plotTheory )
  fig = gcf;
  fig.WindowStyle = dockStyle;
  fprintf('Finished plotBandPhaseSub\n')
  if saveMe
    saveTitle = 'phase_diagram_peak_width';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end
%%
if plotMyIN
  saveExt = 'png';
  xLabel = 'Concentration $$c^*$$';
  yLabel = 'Global nematic order $$ N $$';
  plotOrderParam( bandTable, 'c', 'nAve', xLabel, yLabel )
  fig = gcf;
  fig.WindowStyle = dockStyle;
  fprintf('Finished plot\n')
  if saveMe
    saveTitle = 'in_trans';
    saveas( gcf, saveTitle, saveExt )
    fprintf('Saved\n')
  end
end
%%
