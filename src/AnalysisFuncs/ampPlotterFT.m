function ampPlotterFT(FTmat2plot, FTind2plot, TimeRec, kx0, ky0, km0)

% keyboard

figure()

for i = 1:8
    subplot(2,4,i)
    [Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(i,:) ) ,TimeRec, imag( FTmat2plot(i,:) ) );
    titstr = sprintf('(%d, %d, %d)', ...
        FTind2plot(i,1) - kx0, FTind2plot(i,2) - ky0, FTind2plot(i,3) - km0 );
    title(titstr)
    xlabel('Time')
    
    if(i == 1 || i == 5 ); ylabel(Ax(1),' real Amp'); end
    if(i == 4 || i == 8 ); ylabel(Ax(2),' imag Amp'); end
    
end


% keyboard
end
