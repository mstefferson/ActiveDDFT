function bandPlotFWHM( p, output )
% fwhm
figure()
% C fwhm
subplot(3,1,1);
ylab = 'C';
plot( p.val, output.cFDWHD );
xlabel(p.name); ylabel( [ylab ' FWHM']);
title( [ ylab ' vs ' p.name ] );
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.cFDWHD( output.steady == 0), 'r')

% P fwhm
subplot(3,1,2);
ylab = 'P';
plot( p.val, output.pFDWHD );
xlabel(p.name); ylabel( [ylab ' FWHM']);
title( [ ylab ' vs ' p.name ] );
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.pFDWHD( output.steady == 0), 'r')
% N fwhm
subplot(3,1,3);
ylab = 'N';
plot( p.val, output.nFDWHD );
xlabel(p.name); ylabel( [ylab ' FWHM']);
title( [ ylab ' vs ' p.name ] );
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.nFDWHD( output.steady == 0), 'r')

