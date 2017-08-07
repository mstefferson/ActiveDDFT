% Plot max scaled concentration vs time
function plotMaxCvsTime( C_rec, b, time, saveTag )
if nargin == 3
  saveFlag = 0;
else
  saveFlag = 1;
end
% Find max vs time
Nt = length(time);
% reshape
Cmax = reshape( max( max( C_rec )  ), [1 Nt] );
Cmin = reshape( min( min( C_rec )  ), [1 Nt] );
% set-up fiugure and plot
figure()
%C
Cmin2PLot = Cmin .* b;
Cmax2PLot = Cmax .* b;
[ax] = plotyy( time, Cmax2PLot, time, Cmin2PLot  );
xlabel('time'); ylabel( ax(1),'C max'); ylabel( ax(2),'C min');
ax(1).YLim = [ min(Cmin2PLot) max(Cmax2PLot) ];
ax(2).YLim = [ min(Cmin2PLot) max(Cmax2PLot) ];
fixYLimits( Cmax2PLot, ax(1) )
fixYLimits( Cmin2PLot, ax(2) )
title('Max scaled C');
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
end
