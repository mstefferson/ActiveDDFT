function plotMaxCvsTime( C_rec, b, time, saveTag )
if nargin == 3
  saveFlag = 0;
else
  saveFlag = 1;
end
% Find max vs time
Nt = length(time);
% reshape
C = reshape( max( max( C_rec )  ), [1 Nt] );
% set-up fiugure and plot
figure()
%C
plot( time, C ./ b );
xlabel('time'); ylabel('C'); title('Max scaled C');
% Save it
if saveFlag
  savefig( gcf, [saveTag '.fig'] );
  saveas( gcf, [saveTag '.jpg'], 'jpg')
end
end
