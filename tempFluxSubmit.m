%%
rho = Den_rec(:,:,:,end);
rho_FT = DenFT_rec(:,:,:,end);
C = C_rec(:,:,end);

%%
fluxCmpr( rho, rho_FT, C, flags, systemObj, particleObj, gridObj );