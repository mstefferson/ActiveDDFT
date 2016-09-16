function plotMaxOPvsTime( C_rec, P_rec, N_rec, time, saveTag )

if nargin == 4
  saveFlag = 0;
else
  saveFlag = 1;
end

% Find max vs time
Nt = length(time);

C = reshape( max( max( C_rec )  ), [1 Nt] );
P = reshape( max( max( P_rec ) ), [1 Nt] );
N = reshape( max( max( N_rec ) ), [1 Nt] );


figure()

subplot(1,3,1);
plot( time, C ./ pi );
xlabel('time'); ylabel('C'); title('Max scaled C');


subplot(1,3,2);
plot( time, P );
xlabel('time'); ylabel('Max P'); title('Max Polar Order');

subplot(1,3,3);
plot( time, P );
xlabel('time'); ylabel('Max N'); title('Max Nematic Order');

if saveFlag
  savefig( gcf, [saveTag '.fig'] );
  saveas( gcf, [saveTag '.jpg'], 'jpg')
end
