es1 = 1;
ls1 = 1;
n1 = 1280;
l1 = 20;
n2 = 1280;
l2 = 20;
c = 1;

dk = 2*pi / l1;
kmax = pi*n1/l1;
k = dk*( -n1/2:n1/2-1 );
paramVec = [n1 n2 l1 l2 es1 ls1 c];

[v, vFt] = tempGaussian( es1, ls1, n1, l1, n2, l2 );

%%
vFt = real(vFt);
imagesc( vFt )
plot( k .* ls1, vFt(n1/2+1,:) )

 