function ampPlotterFT(FTmat2plot, TimeRec, Nx, Ny, Nm, bc, vD,SaveMe,trial)

ParamStrNx = sprintf('Nx = %d', Nx);
ParamStrNy = sprintf('Ny = %d', Ny);
ParamStrNm = sprintf('Nm = %d', Nm);
ParamStrBc = sprintf('Nm = %.2f', bc);
ParamStrVd = sprintf('Nm = %.2f', vD);

figure()

subplot(2,4,1)

% keyboard
[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(1,:) ) ,TimeRec, imag( FTmat2plot(1,:) ) );
title('Mode (0,0,1)')
xlabel('Time')
ylabel(Ax(1),' real Amp')
textbp(ParamStrNx)

subplot(2,4,2)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(2,:) ) ,TimeRec, imag( FTmat2plot(3,:) ) );
title('Mode (1,0,1)')
xlabel('Time')
ylabel(Ax(1),' real Amp')
textbp(ParamStrNy)

subplot(2,4,3)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(3,:) ) ,TimeRec, imag( FTmat2plot(3,:) ) );
title('Mode (0,1,1)')
xlabel('Time')
ylabel(Ax(1),' real Amp')
textbp(ParamStrNm)

subplot(2,4,4)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(4,:) ) ,TimeRec, imag( FTmat2plot(4,:) ) );
title('Mode (1,1,1)')
xlabel('Time')
ylabel(Ax(1),' real Amp')


subplot(2,4,5)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(5,:) ) ,TimeRec, imag( FTmat2plot(5,:) ) );
title('Mode (0,0,2)')
xlabel('Time')
ylabel(Ax(1),' real Amp')
textbp(ParamStrVd)

subplot(2,4,6)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(6,:) ) ,TimeRec, imag( FTmat2plot(6,:) ) );
title('Mode (1,0,2)')
xlabel('Time')
ylabel(Ax(1),' real Amp')
textbp(ParamStrBc)

subplot(2,4,7)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(7,:) ) ,TimeRec, imag( FTmat2plot(7,:) ) );
title('Mode (0,1,2)')
xlabel('Time')
ylabel(Ax(1),' real Amp')

subplot(2,4,8)

[Ax, ~, ~] = plotyy( TimeRec, real( FTmat2plot(8,:) ) ,TimeRec, imag( FTmat2plot(8,:) ) );
title('Mode (1,1,2)')
xlabel('Time')
ylabel(Ax(1),' real Amp')

if SaveMe
figtl = sprintf('AmpFT%d',trial);
savefig(gcf,figtl)
saveas(gcf, figtl,'jpg')
end