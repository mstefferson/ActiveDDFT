%% 
cWant = [1.55 1.60 1.65 1.70 1.75 1.80 1.875];
fWant = 10;
% cWant = [1.60];
% fWant = [0 5 7.5 10 15];
%%
figure()
hold on
legCell = cell( 1, length( cWant ) * length( fWant ) );
counter = 1;
for ii = 1:length(cWant)
  for jj = 1:length(fWant)
    cWantTemp = cWant(ii);
    fWantTemp = fWant(jj);
    logInds = bandTable.c == cWantTemp & bandTable.fd == fWantTemp;
    cTemp = bandTable.cSlice{ logInds };
    cMax = bandTable.cMax( logInds );
    cFWHM = bandTable.cFWHM( logInds );
    fprintf( 'c = %f fd = %f cMax = %f, cFWHM = %f\n', ...
      cWantTemp, fWantTemp, cMax ./ cWantTemp, cFWHM )
    legCell{counter} = [ 'c=' num2str(cWantTemp,'%.3f') ...
      ' fd=' num2str(fWantTemp,'%.1f') ];
    plot( cTemp )
    counter = counter + 1;
  end
end
legend( legCell, 'location', 'best' )
hold off