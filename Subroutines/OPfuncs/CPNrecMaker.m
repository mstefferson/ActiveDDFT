function [OrderParamObj] = ...
    CPNrecMaker(Nx,Ny,TimeRecVec,GridObj,Density_rec)

% Number of frames and time vec
OrderParamObj.nFrames    = length(TimeRecVec) ; 
OrderParamObj.TimeRecVec = TimeRecVec ; 
 % Concentration
OrderParamObj.C_rec      = zeros(Nx,Ny,nFrames);
% Polar Order Parameter
OrderParamObj.POP_rec    = zeros(Nx,Ny,nFrames);
% 1st moment orientation of distribution-x,y dir
OrderParamObj.POPx_rec = zeros(Nx,Ny,nFrames);
OrderParamObj.POPy_rec = zeros(Nx,Ny,nFrames);
% Nematic order parameter-(max eigenvalue of NOP)
OrderParamObj.NOP_rec    = zeros(Nx,Ny,nFrames);
% Nematic alignment director-x,y dir (eigenvector of nematic order parameter)
OrderParamObj.NADx_rec   = zeros(Nx,Ny,nFrames);
OrderParamObj.NADy_rec   = zeros(Nx,Ny,nFrames);


for ii = 1:nFrames
    [C,POP,POPx,POPy,NOP,NADx,NADy] = ...
        OpCPNCalc(Nx,Ny,Density_rec(:,:,:,ii),...
        GridObj.phi,GridObj.x,GridObj.y,GridObj.phi3D);
    
    OrderParamObj.C_rec(:,:,ii) = real(C);
    
    OrderParamObj.POP_rec(:,:,ii)    = real(POP);
    OrderParamObj.POPx_rec(:,:,ii) = real(POPx);
    OrderParamObj.POPy_rec(:,:,ii) = real(POPy);
    
    OrderParamObj.NOP_rec(:,:,ii)  = real(NOP);
    OrderParamObj.NADx_rec(:,:,ii) = real(NADx);
    OrderParamObj.NADy_rec(:,:,ii) = real(NADy);
%     keyboard 
end

end
