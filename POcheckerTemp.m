%%
Delta = 17;
Tind  = 129-Delta;
Mind  = 65;
Bind  = 1 + Delta;

rho = reshape( DenRecObj.Density_rec(:,:,:,end), [128,128,128] );

C = reshape( OrderParamObj.C_rec(:,:,end), [128,128] );
PO = reshape( OrderParamObj.POP_rec(:,:,end), [128,128] );
NO = reshape( OrderParamObj.NOP_rec(:,:,end), [128,128] );

POx = reshape( OrderParamObj.nx_POP_rec(:,:,end), [128,128] );
POy = reshape( OrderParamObj.ny_POP_rec(:,:,end), [128,128] );

NOx = reshape( OrderParamObj.NADx_rec(:,:,end), [128,128] );
NOy = reshape( OrderParamObj.NADy_rec(:,:,end), [128,128] );


phi = GridObj.phi;
Ny = trapz_periodic(GridObj.phi, sin(GridObj.phi3D) .* rho,3);
Nx = trapz_periodic(GridObj.phi, cos(GridObj.phi3D) .* rho,3);

% Surface Polar Order Directors
figure()
subplot(2,2,1)
pcolor(Nx'); colorbar;
shading('interp');
xlabel('x');xlabel('y');
title('Nx');

subplot(2,2,2)
pcolor(Ny');colorbar;
shading('interp');
xlabel('x');xlabel('y');
title('Ny');

% Surface Polar Order Directors
subplot(2,2,3)
pcolor(C'); colorbar;
shading('interp');
xlabel('x');xlabel('y');
title('C');

subplot(2,2,4)
pcolor(NO');colorbar;
shading('interp');
xlabel('x');xlabel('y');
title('NO');


% Distributions
figure()
DistT = reshape( rho(1,Tind,:), [1, 128]  );
DistM = reshape( rho(1,Mind,:), [1, 128]  );
DistB = reshape( rho(1,Bind,:), [1, 128]  );

plot( GridObj.phi, DistT, GridObj.phi, DistM,GridObj.phi, DistB)
xlabel('\phi'); ylabel('f');
legend('top','middle','bottom')

nxT = trapz_periodic(phi, cos(phi).* DistT);  %Polar order parameter in x-direction
nyT = trapz_periodic(phi, sin(phi).*DistT);  %Polar order parameter in y-direction

nxTnorm = nxT ./  C(1,Tind);
nyTnorm = nyT ./  C(1,Tind);

nxM = trapz_periodic(phi, cos(phi).* DistM); %Polar order parameter in x-direction
nyM = trapz_periodic(phi, sin(phi).*DistM);  %Polar order parameter in y-direction

nxMnorm = nxM ./ C(1,Mind);
nyMnorm = nyM ./ C(1,Mind);

nxB = trapz_periodic(phi, cos(phi).* DistB);  %Polar order parameter in x-direction
nyB = trapz_periodic(phi, sin(phi).*DistB); %Polar order parameter in y-direction

nxBnorm = nxB ./ C(1,Bind);
nyBnorm = nyB ./ C(1,Bind);

figure()

subplot(2,2,1)
pcolor(C'); shading('interp')
xlabel('x');xlabel('y');
title('C');

subplot(2,2,2)
pcolor(PO'); shading('interp')
xlabel('x');xlabel('y');
title('PO');

subplot(2,2,3)
plot( GridObj.phi, DistT, GridObj.phi, DistM,GridObj.phi, DistB)
legend('top','middle','bottom','location','best')
xlabel('\phi');xlabel('f');
title('Dist');

subplot(2,2,4)
pcolor(NO'); shading('interp')
xlabel('x');xlabel('y');
title('NO');

figure()
subplot(1,3,1)
plot(GridObj.y, C(1,:) ./ max(max( C ) ),...
    GridObj.y, PO(1,:) ./ max(max( PO ) ), ...
    GridObj.y, NO(1,:) ./ max(max( NO ) ) );
xlabel('y');xlabel('Normalized value');
title('OPs')
legend('C','PO','NO','location','best')

subplot(1,3,2)
plot(GridObj.y, POy(1,:))
title('PO y dir')
xlabel('y');xlabel('POy');


subplot(1,3,3)
plot(GridObj.y, POy(1,:) ./ max( POy(1,: ) ), ...
    GridObj.y, C(1,:) ./ max( C(1,: ) ),...
    GridObj.y, NO(1,:) ./ max( NO(1,: ) ))
xlabel('y');xlabel('Normalized value');
title('C, NO, PO y dir')
legend('POy','C','NO','location','best')



