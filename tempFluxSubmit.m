%%
rho = denRecObj.rhoFinal;
rhoFT = fftshift(fftn( rho ) );
C = C_rec(:,:,end);

%%
fluxCmpr( rho, rhoFT, C, flags, systemObj, particleObj, gridObj );