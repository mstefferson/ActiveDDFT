function [OrderParamObj] = ...
    CPNrecMaker(n1,n2,TimeRecVec,Density_rec,...
    phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d )

% Number of frames and time vec
nFrames = length(TimeRecVec);
OrderParamObj.nFrames    = nFrames; 
OrderParamObj.TimeRecVec = TimeRecVec ; 
 % Concentration
OrderParamObj.C_rec      = zeros(n1,n2,nFrames);
% Polar Order Parameter
OrderParamObj.POP_rec    = zeros(n1,n2,nFrames);
% 1st moment orientation of distribution-x,y dir
OrderParamObj.POPx_rec = zeros(n1,n2,nFrames);
OrderParamObj.POPy_rec = zeros(n1,n2,nFrames);
% Nematic order parameter-(max eigenvalue of NOP)
OrderParamObj.NOP_rec    = zeros(n1,n2,nFrames);
% Nematic alignment director-x,y dir (eigenvector of nematic order parameter)
OrderParamObj.NOPx_rec   = zeros(n1,n2,nFrames);
OrderParamObj.NOPy_rec   = zeros(n1,n2,nFrames);


for ii = 1:nFrames
    [C,POP,POPx,POPy,NOP,NOPx,NOPy] = ...
        OpCPNCalc(n1,n2,Density_rec(:,:,:,ii),...
        phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
    
    OrderParamObj.C_rec(:,:,ii) = real(C);
    
    OrderParamObj.POP_rec(:,:,ii)    = real(POP);
    OrderParamObj.POPx_rec(:,:,ii) = real(POPx);
    OrderParamObj.POPy_rec(:,:,ii) = real(POPy);
    
    OrderParamObj.NOP_rec(:,:,ii)  = real(NOP);
    OrderParamObj.NOPx_rec(:,:,ii) = real(NOPx);
    OrderParamObj.NOPy_rec(:,:,ii) = real(NOPy);
%     keyboard 
end

end
