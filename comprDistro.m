function comprDistro( nVec )
% list of files to examine
bcVec = [1.35 1.45 1.55 1.65];
% file Ids
fileStartId = 'Hr_rods_mayer_diag1_N';
fileEnd = { ...
  '_ls1010_bc1.35_vD0_IC1_SM6_t01.02/',...
  '_ls1010_bc1.45_vD0_IC1_SM6_t01.02/',...
  '_ls1010_bc1.55_vD0_IC1_SM6_t01.02/',...
  '_ls1010_bc1.65_vD0_IC1_SM6_t01.02/',...
  '_ls1010_bc1.35_vD0_IC1_SM6_t01.01/',...
  '_ls1010_bc1.45_vD0_IC1_SM6_t01.01/',...
  '_ls1010_bc1.55_vD0_IC1_SM6_t01.01/',...
  '_ls1010_bc1.65_vD0_IC1_SM6_t01.01/'};
% file path
filePath = './analyzedfiles/analyzeMe/testIN/';
% call functions
comprDistroBcPanel( filePath, fileStartId, fileEnd, nVec, bcVec );
comprDistroNPanel( filePath, fileStartId, fileEnd, nVec, bcVec );
end
%
% functions
% comprDistroIsoNem( filePath, fileStartId, fileEnd, nVec )
function [fig1] = comprDistroBcPanel( filePath, fileStartId, fileEnd, nVec, bcVec )
numNs = length( nVec );
numBcs = length( bcVec );
% set up fig
fig1 = figure();
axVec = zeros( 2*numBcs, 1 );
for ii = 1:numBcs
  axVec( ii ) = subplot(2, numBcs,ii);
  hold
  xlabel('$$ \phi $$');
  ylabel('$$ f $$');
  title( ['old: bc = ' num2str( bcVec(ii) )] );
  if bcVec(ii) < 1.50; ax=gca; ax.YLim=[0 1.1]; end
  axVec( numBcs + ii ) = subplot(2,numBcs, numBcs + ii);
  hold
  xlabel('$$ \phi $$');
  ylabel('$$ f $$');
  title( ['new: bc = ' num2str( bcVec(ii) )] );
  if bcVec(ii) < 1.50; ax=gca; ax.YLim=[0 1.1]; end
end
% legend
legcell = cell( numNs, 1 );
% loop and plot
for ii = 1:numNs
  nTemp = nVec(ii);
  xLim = 2*pi /  nVec(ii) * [0 nTemp - 1];
  insertStr = sprintf('%d%d%d', nTemp,nTemp,nTemp);
  for jj = 1:numBcs
    % load it
    fileId = [fileStartId insertStr fileEnd{jj}];
    fullPath = [filePath fileId];
    % get distro
    [f, phi] = getDist( fullPath,  nTemp );
    % plot it
    subplot(2,numBcs,jj);
    ax = gca;
    ax.XLim = xLim;
    plot(phi, f)
    % load it
    fileId = [fileStartId insertStr fileEnd{jj+numBcs}];
    fullPath = [filePath fileId];
    % get distro
    [f, phi] = getDist( fullPath, nVec(ii) );
    % plot it
    subplot(2,numBcs,jj+numBcs);
    ax = gca;
    ax.XLim = xLim;
    plot(phi, f)
  end
  legcell{ii} = ['N=' num2str( nVec(ii) )];
end

% put legend on
leg = legend( axVec(end), legcell, 'location', 'best');
leg.FontSize = 8;
end
%
% comprDistroAllBc( filePath, fileStartId, fileEnd, nVec, bcVec )
function [fig1] = comprDistroNPanel( filePath, fileStartId, fileEnd, nVec, bcVec )
% get lengths
numNs = length( nVec );
numBcs = length( bcVec );
% get figs locked and loaded
fig1 = figure();
axVec = zeros( 2*numNs, 1 );
for ii = 1:numNs
  nTemp = nVec(ii);
  xLim = 2*pi /  nVec(ii) * [0 nTemp - 1];
  axVec( ii ) = subplot(2, numNs,ii);
  ax = gca;
  ax.XLim = xLim;
  hold
  xlabel('$$ \phi $$');
  ylabel('$$ f $$');
  title( ['old: N = ' num2str( nTemp )] );
  axVec( numNs + ii ) = subplot(2,numNs, numNs + ii);
  ax = gca;
  ax.XLim = xLim;
  hold
  xlabel('$$ \phi $$');
  ylabel('$$ f $$');
  title( ['new: N = ' num2str( nTemp )] );
end
% legend
legcell = cell( numBcs, 1 );
% loop and plot
for ii = 1:numNs
  nTemp = nVec(ii);
  insertStr = sprintf('%d%d%d',nTemp,nTemp,nTemp); 
  for jj = 1:numBcs
    fileId = [fileStartId insertStr fileEnd{jj}];
    fullPath = [filePath fileId];
    % get distro
    [f, phi] = getDist( fullPath, nTemp );
    % plot it
    subplot(2,numNs,ii);
    plot(phi, f)
    fileId = [fileStartId insertStr fileEnd{jj+numBcs}];
    fullPath = [filePath fileId];
    % get distro
    [f, phi] = getDist( fullPath, nTemp);
    % plot it
    subplot(2,numNs,ii+numNs);
    ax = gca;
    plot(phi, f)
    ax.XLim = xLim;
    if ii == 1
      legcell{jj} = ['c*=' num2str( bcVec(jj) )];
    end
  end
end
% put legend on
leg = legend( axVec(end), legcell, 'location', 'best');
leg.FontSize = 8;
end
% getDist
function [f,phi] = getDist( fullPath, n )
% load params and check for rhoFinal
paramFile = dir( [fullPath 'params*' ]);
if exist([fullPath paramFile.name],'file' )
  load( [fullPath paramFile.name] );
  if isfield( denRecObj, 'rhoFinal' )
    rhoFinal = denRecObj.rhoFinal;
  else
    rhoFinalFile = dir( [fullPath 'rhoFinal*'] );
    load( [fullPath rhoFinalFile.name]);
    rhoFinal = rho;
  end
  % get distrbution
  f = reshape( rhoFinal(1,1,:), [1 n] );
  % get phi
  phi = 2*pi / n * (0:n-1);
else
  f = []; phi=[];
end
end
