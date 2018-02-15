d0 = 1;
rhoScale = 1;
c_max = 6;
conc = linspace( 0, c_max );


% calcDnl with rods function.
dnl_rods = d0 + d0 .* (  1 ./ ( 1  +  conc / rhoScale  ) .^2  - 1 );

dnl_rods_wrong = d0 + d0 .* (  1 ./ ( 1  + ( conc / rhoScale  ) .^2 ) - 1 );

% calcDnl with exp function.
dnl_exp = d0 + d0 .* (  exp( -conc / rhoScale ) - 1 );

% calcDnl with tanh function.
dnl_tanh = d0 -d0 .* tanh( conc / rhoScale );

figure()
plot( conc, dnl_tanh, conc, dnl_exp, conc, dnl_rods )
xlabel('c')
ylabel('D')
legend('tanh', 'exp', 'rods', 'location', 'best')

figure()
plot( conc, dnl_rods, conc, dnl_rods_wrong )
xlabel('c')
ylabel('D')

figure()
slopebig = 10000;
cin = 1.5;
p = plot( conc, dnl_rods );
p.LineStyle = '-';
hold on
p = plot( conc, dnl_rods_wrong );
p.LineStyle = '-';
p = plot(conc, slopebig*(conc-cin) );
p.LineStyle = ':';
hold off
ax = gca;
ax.YLim = [0 1];
xlabel('$$ c $$')
ylabel('$$ D $$')
l = legend('$$ D $$', '$$ D_{wrong} $$', '$$ c_{IN} $$');
l.Interpreter = 'latex';

figure()
slopebig = 10000;
cin = 1.5;
p = plot( conc, dnl_rods );
p.LineStyle = '-';
hold on
p = plot(conc, slopebig*(conc-cin) );
p.LineStyle = ':';
hold off
ax = gca;
ax.YLim = [0 1];
xlabel('$$ c $$')
ylabel('$$ D $$')
l = legend('$$ D $$', '$$ c_{IN} $$');
l.Interpreter = 'latex';
axis square