function [OrderParamObj] = ...
    CPNrecMaker(Nx,Ny,TimeRecVec,GridObj,Density_rec,feq)

nFrames    = length(TimeRecVec)  ;                 % Number of frames
C_rec      = zeros(Nx,Ny,nFrames);                 % Concentration
POP_rec    = zeros(Nx,Ny,nFrames);                 % Polar Order Parameter
nx_POP_rec = zeros(Nx,Ny,nFrames);                 % 1st moment orientation of distribution-x dir
ny_POP_rec = zeros(Nx,Ny,nFrames);                 % 1st moment orientation of distribution-x dir
NOP_rec    = zeros(Nx,Ny,nFrames);                 % Nematic order parameter-(max eigenvalue of NOP)
NADx_rec   = zeros(Nx,Ny,nFrames);                 % Nematic alignment director-x dir (eigenvector of nematic order parameter)
NADy_rec   = zeros(Nx,Ny,nFrames);                 % Nematic alignment director-y dir (eigenvector of nematic order parameter)


for ii = 1:nFrames
    [C,POP,nx_POP,ny_POP,NOP,NADx,NADy] = ...
        OpCPNCalc(Nx,Ny,Density_rec(:,:,:,ii),...
        GridObj.phi,GridObj.x,GridObj.y,GridObj.phi3D);
    
    C_rec(:,:,ii) = real(C);
    
    POP_rec(:,:,ii)    = real(POP);
    nx_POP_rec(:,:,ii) = real(nx_POP);
    ny_POP_rec(:,:,ii) = real(ny_POP);
    
    NOP_rec(:,:,ii)  = real(NOP);
    NADx_rec(:,:,ii) = real(NADx);
    NADy_rec(:,:,ii) = real(NADy);
%     keyboard 
end

% keyboard
[~,~,~,~,NOPeq,~,~] = ...
        OpCPNCalc(1,1,feq,...
        GridObj.phi,1,1,GridObj.phi3D);
    
OrderParamObj = struct('nFrames',nFrames,'TimeRec',TimeRecVec,...
    'C_rec',C_rec, 'POP_rec',POP_rec,...
    'nx_POP_rec',nx_POP_rec,'ny_POP_rec',ny_POP_rec,...
                       'NOP_rec',NOP_rec,'NADx_rec',...
                       NADx_rec,'NADy_rec',NADy_rec,...
                       'NOPeq',NOPeq);

end