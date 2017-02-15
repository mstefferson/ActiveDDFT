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
% set-up fiugure and plot
figure()
%C
subplot(1,3,1);
plot( time, C ./ b );
xlabel('time'); ylabel('C'); title('Max scaled C');
%P
subplot(1,3,2);
plot( time, P );
xlabel('time'); ylabel('Max P'); title('Max Polar Order');
%N
subplot(1,3,3);
plot( time, N );
xlabel('time'); ylabel('Max N'); title('Max Nematic Order');
% Save it
if saveFlag
  savefig( gcf, [saveTag '.fig'] );
  saveas( gcf, [saveTag '.jpg'], 'jpg')
end
