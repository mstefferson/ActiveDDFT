function ModesVsTimePlotterKm(ampl_record,k2plotInd,TimeRecVec, n3odes,kxholder,kyholder,...
    n1,n2,n3)
    
    figure(9)
    % Plot the mode amplitudes except for the k =
    NumTicks = 3; % Number of Tick marks on the y axis
    
    for j = 1:(n3odes + 1 ) / 2
        % Make a vector of single mode amplitude throughout time.
        %Make them all positive so they all look similiar. Sign doesn't matter anyway
        
        if k2plotInd(j) == (n3/2 + 1)
            
            subplot( (n3odes+1) / 2, 2,[1,2] )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plotInd(j),:) ),...
                TimeRecVec,imag( ampl_record(k2plotInd(j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                kxholder - ( n1 / 2 + 1 ), ...
                kyholder - ( n2 / 2 + 1 ), ...
                k2plotInd(j)- ( n3 / 2 + 1 ) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            xlabel(haxes(2),'time')
            set(haxes(1),'YTick', mean( real( ampl_record(k2plotInd(j),:) ) ) )
            set(haxes(2),'YTick',mean( imag( ampl_record(k2plotInd(j),:) ) ))
            title(titlestr)
            
        else
            
%             keyboard
            %plot negative mode first
            subplot( (n3odes+1) / 2, 2, (n3odes + 3) -  2*j )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plotInd(j),:) ),...
                TimeRecVec,imag( ampl_record(k2plotInd(j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                kxholder - ( n1 / 2 + 1 ), ...
                kyholder - ( n2 / 2 + 1 ), ...
                k2plotInd(j)- ( n3 / 2 + 1 ) );
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
            subplot( (n3odes+1) / 2, 2, (n3odes + 2) -  2*j )
            [haxes,hline1,hline2] = plotyy(TimeRecVec,real( ampl_record(k2plotInd(n3odes + 1 - j),:) ),...
                TimeRecVec,imag( ampl_record(k2plotInd(n3odes + 1 - j),:) ) );
            
            titlestr = sprintf('kx = %i ky = %i km = %i',...
                kxholder - ( n1 / 2 + 1 ), ...
                kyholder - ( n2 / 2 + 1 ), ...
                k2plotInd(n3odes + 1 - j) - (n3/2 + 1) );
            ylabel(haxes(1),'Re \{F_k\}')
            ylabel(haxes(2),'Im \{F_k\}')
            %             xlabel(haxes(2),'time')
            
            
            WantedYTicks = linspace(...
                min( real( ampl_record(k2plotInd(n3odes + 1 - j),:) ) ) , ...
                max( real( ampl_record(k2plotInd(n3odes + 1 - j),:) ) ),  NumTicks);
            set(haxes(1),'YTick', WantedYTicks)
            
            WantedYTicks = linspace(...
                min( imag( ampl_record(k2plotInd(n3odes + 1 - j),:) ) ) , ...
                max( imag( ampl_record(k2plotInd(n3odes + 1 - j),:) ) ),  NumTicks);
            set(haxes(2),'YTick',WantedYTicks)
            
            title(titlestr)
        end %end if k =0 statement
    end %end mode loop
    
end % end AllKsVsTime
