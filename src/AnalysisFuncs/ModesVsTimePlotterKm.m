function ModesVsTimePlotterKm(ampl_record,k2plotInd,TimeRecVec, Nmodes,kxholder,kyholder,...
    Nx,Ny,Nm)
    
    figure(9)
    % Plot the mode amplitudes except for the k =
    NumTicks = 3; % Number of Tick marks on the y axis
    
    for j = 1:(Nmodes + 1 ) / 2
        % Make a vector of single mode amplitude throughout time.
        %Make them all positive so they all look similiar. Sign doesn't matter anyway
        
        if k2plotInd(j) == (Nm/2 + 1)
            
            subplot( (Nmodes+1) / 2, 2,[1,2] )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plotInd(j),:) ),...
                TimeRecVec,imag( ampl_record(k2plotInd(j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                kxholder - ( Nx / 2 + 1 ), ...
                kyholder - ( Ny / 2 + 1 ), ...
                k2plotInd(j)- ( Nm / 2 + 1 ) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            xlabel(haxes(2),'time')
            set(haxes(1),'YTick', mean( real( ampl_record(k2plotInd(j),:) ) ) )
            set(haxes(2),'YTick',mean( imag( ampl_record(k2plotInd(j),:) ) ))
            title(titlestr)
            
        else
            
%             keyboard
            %plot negative mode first
            subplot( (Nmodes+1) / 2, 2, (Nmodes + 3) -  2*j )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plotInd(j),:) ),...
                TimeRecVec,imag( ampl_record(k2plotInd(j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                kxholder - ( Nx / 2 + 1 ), ...
                kyholder - ( Ny / 2 + 1 ), ...
                k2plotInd(j)- ( Nm / 2 + 1 ) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            %             xlabel(haxes(2),'time')
            WantedYTicks = linspace(...
                min( real( ampl_record(k2plotInd(j),:) ) ) , ...
                max( real( ampl_record(k2plotInd(j),:) ) ),  NumTicks);
            set(haxes(1),'YTick', WantedYTicks)
            
            WantedYTicks = linspace(...
                min( imag( ampl_record(k2plotInd(j),:) ) ) , ...
                max( imag( ampl_record(k2plotInd(j),:) ) ),  NumTicks);
            set(haxes(2),'YTick',WantedYTicks)
            title(titlestr)
            
            %Then plus
            subplot( (Nmodes+1) / 2, 2, (Nmodes + 2) -  2*j )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plotInd(Nmodes + 1 - j),:) ),...
                TimeRecVec,imag( ampl_record(k2plotInd(Nmodes + 1 - j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                kxholder - ( Nx / 2 + 1 ), ...
                kyholder - ( Ny / 2 + 1 ), ...
                k2plotInd(Nmodes + 1 - j) - (Nm/2 + 1) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            %             xlabel(haxes(2),'time')
            
            
            WantedYTicks = linspace(...
                min( real( ampl_record(k2plotInd(Nmodes + 1 - j),:) ) ) , ...
                max( real( ampl_record(k2plotInd(Nmodes + 1 - j),:) ) ),  NumTicks);
            set(haxes(1),'YTick', WantedYTicks)
            
            WantedYTicks = linspace(...
                min( imag( ampl_record(k2plotInd(Nmodes + 1 - j),:) ) ) , ...
                max( imag( ampl_record(k2plotInd(Nmodes + 1 - j),:) ) ),  NumTicks);
            set(haxes(2),'YTick',WantedYTicks)
            
            title(titlestr)
        end %end if k =0 statement
    end %end mode loop
    
end % end AllKsVsTime
