function ampPlotterFT(FTmat2plot, FTind2plot, TimeRec, Nx, Ny, Nm, bc, vD, SaveMe,trial)

kx0 = Nx / 2 + 1;
ky0 = Ny / 2 + 1;
km0 = Nm / 2 + 1;

% keyboard

ParamStrNx = sprintf('Nx:%d', Nx);
ParamStrNy = sprintf('Ny:%d', Ny);
ParamStrNm = sprintf('Nm:%d', Nm);
ParamStrBc = sprintf('bc:%.2f', bc);
ParamStrVd = sprintf('vD:%.2f', vD);
ParamStrTr = sprintf('t:%d', trial);

ParamStrCell = {ParamStrNx;  ParamStrNy; ParamStrNm;' ';...
                         ParamStrTr;ParamStrBc; ParamStrVd;' ' };
    

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
    
%     keyboard
%     textbp(ParamStrCell{i})
%     textbp('helloasdfadfa')
end

if SaveMe
figtl = sprintf('AmpFT%d',trial);
savefig(gcf,figtl)
saveas(gcf, figtl,'jpg')
end

% keyboard
end