function ampPlotterFT(FTmat2plot, FTind2plot, TimeRec, Nx, Ny, Nm, bc, vD, SaveMe,trial)

kx0 = Nx / 2 + 1;
ky0 = Ny / 2 + 1;
km0 = Nm / 2 + 1;

% keyboard

ParamStrNx = sprintf('Nx = %d', Nx);
ParamStrNy = sprintf('Ny = %d', Ny);
ParamStrNm = sprintf('Nm = %d', Nm);
ParamStrBc = sprintf('bc = %.2f', bc);
ParamStrVd = sprintf('vD = %.2f', vD);

ParamStrCell = {ParamStrNx;  ParamStrNy; ParamStrNm;' ';...
                         ParamStrBc; ParamStrVd;' ';' ' };
    

figure()

% keyboard
for i = 1:8
    subplot(2,4,i)
    [Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(i,:) ) ,TimeRec, imag( FTmat2plot(i,:) ) );
    titstr = sprintf('Mode (%d, %d, %d)', ...
        FTind2plot(i,1) - kx0, FTind2plot(i,2) - ky0, FTind2plot(i,3) - km0 );
    title(titstr)
    xlabel('Time')
    
    if(i == 1 || i == 5 ); ylabel(Ax(1),' real Amp'); end
    if(i == 4 || i == 8 ); ylabel(Ax(2),' imag Amp'); end
    
    textbp(ParamStrCell{i})
end
textbp(ParamStrNy)

textbp(ParamStrNm)

subplot(2,4,5)

textbp(ParamStrVd)

textbp(ParamStrBc)

if SaveMe
figtl = sprintf('AmpFT%d',trial);
savefig(gcf,figtl)
saveas(gcf, figtl,'jpg')
end

% keyboard
end