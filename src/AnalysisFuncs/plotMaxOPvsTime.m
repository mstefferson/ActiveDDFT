function plotMaxOPvsTime( C_rec, P_rec, N_rec, b, time, saveTag )
if nargin == 5
  saveFlag = 0;
else
  saveFlag = 1;
end
% Find max vs time
Nt = length(time);
% reshape
C = reshape( max( max( C_rec )  ), [1 Nt] );
P = reshape( max( max( P_rec ) ), [1 Nt] );
N = reshape( max( max( N_rec ) ), [1 Nt] );
% scale C
C = C .* b;
% set-up fiugure and plot
figure()
%% C
ax = subplot(1,3,1);
plot( time, C );
xlabel('time'); ylabel('C'); title('Max scaled C');
fixYLimits( C, ax );
%% P
ax = subplot(1,3,2);
plot( time, P );
xlabel('time'); ylabel('Max P'); title('Max Polar Order');
fixYLimits( P, ax );
%% N
ax = subplot(1,3,3);
plot( time, N );
xlabel('time'); ylabel('Max N'); title('Max Nematic Order');
fixYLimits( N, ax );
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


