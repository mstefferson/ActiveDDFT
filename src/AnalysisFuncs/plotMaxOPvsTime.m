function plotMaxOPvsTime( C_rec, P_rec, N_rec, b, time, saveTag )
if nargin == 5
  saveFlag = 0;
else
  saveFlag = 1;
end
% Find max vs time
Nt = length(time);
% reshape
Cmax = b .* reshape( max( max( C_rec )  ), [1 Nt] );
Pmax = reshape( max( max( P_rec ) ), [1 Nt] );
Nmax = reshape( max( max( N_rec ) ), [1 Nt] );
Cmin = b .* reshape( min( min( C_rec )  ), [1 Nt] );
Pmin = reshape( min( min( P_rec ) ), [1 Nt] );
Nmin = reshape( min( min( N_rec ) ), [1 Nt] );
% set-up fiugure and plot
figure()
%% C
ax = subplot(1,3,1);
plot( time, Cmax );
hold
plot( time, Cmin );
xlabel('time'); ylabel('C'); title('Max scaled C');
% fixYLimits( Cmax, ax );
%% P
ax = subplot(1,3,2);
plot( time, Pmax );
hold
plot( time, Pmin );
xlabel('time'); ylabel('Max P'); title('Max Polar Order');
% fixYLimits( Pmax, ax );
%% N
ax = subplot(1,3,3);
plot( time, Nmax );
hold
plot( time, Nmin );
xlabel('time'); ylabel('Max N'); title('Max Nematic Order');
% fixYLimits( Nmax, ax );
% Save it
if saveFlag
  savefig( gcf, [saveTag '.fig'] );
  saveas( gcf, [saveTag '.jpg'], 'jpg')
end
% subroutines
% fix limits
  function fixYLimits( v, ax )
    scaleEps = 0.01;
    meanVal = mean(v);
    if abs( std( v ) ./ mean(v) ) < scaleEps
      minLim = min(v) - meanVal * scaleEps;
      maxLim = max(v) + meanVal * scaleEps;
      ax.YLim = [ minLim maxLim ];
    end
  end % function
end % max plot


