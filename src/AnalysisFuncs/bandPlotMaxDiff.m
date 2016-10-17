function bandPlotMaxDiff( p, output )
% Max and diff
figure()
% C max and diff
subplot(3,1,1);
ylab = 'C';
[Ax] = plotyy( p.val, output.cMax, p.val, output.cDiff );
xlabel(p.name); ylabel(Ax(1), [ ylab ' Max']); ylabel(Ax(2), [ ylab ' diff']);
legend('Max', 'delta');
title( [ ylab ' vs ' p.name ] );
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.cMax( output.steady == 0), 'r')

% P max and diff
subplot(3,1,2);
ylab = 'P';
[Ax] = plotyy( p.val, output.pMax, p.val, output.pDiff );
xlabel(p.name); ylabel(Ax(1), [ ylab ' Max']); ylabel(Ax(2), [ ylab ' diff']);
title( [ ylab ' vs ' p.name ] );
ylab = 'Max P';
xlabel(p.name); ylabel(ylab);
legend('Max', 'delta');
title( [ ylab ' vs ' p.name ] );
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.pMax( output.steady == 0), 'r')

% N max and diff
subplot(3,1,3);
ylab = 'N';
[Ax] = plotyy( p.val, output.nMax, p.val, output.nDiff );
xlabel(p.name); ylabel(Ax(1), [ ylab ' Max']); ylabel(Ax(2), [ ylab ' diff']);
legend('Max', 'delta');
title( [ ylab ' vs ' p.name ] );
% plot red circles if steady state was no reached
hold on
scatter(p.val( output.steady == 0), output.nMax( output.steady == 0), 'r')

end
