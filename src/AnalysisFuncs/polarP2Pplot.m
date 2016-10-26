function polarP2Pplot( p, output )

figure()
plot(p.val, output.pP2P, 'x' );
title('Distance between polar peaks');
xlabel(p.name); ylabel( ['peak distance']);
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.pP2P( output.steady == 0), 'r')