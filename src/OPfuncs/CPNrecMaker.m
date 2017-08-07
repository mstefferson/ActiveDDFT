function [orderParamObj] = ...
    CPNrecMaker(n1,n2,n3,timeRecVec,Density_rec,...
    phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d )

% Number of frames and time vec
nFrames = length(timeRecVec);
orderParamObj.nFrames    = nFrames; 
orderParamObj.timeRecVec = timeRecVec ; 
 % Concentration
orderParamObj.C_rec      = zeros(n1,n2,nFrames);
% Polar Order Parameter
orderParamObj.POP_rec    = zeros(n1,n2,nFrames);
% 1st moment orientation of distribution-x,y dir
orderParamObj.POPx_rec = zeros(n1,n2,nFrames);
orderParamObj.POPy_rec = zeros(n1,n2,nFrames);
% Nematic order parameter-(max eigenvalue of NOP)
orderParamObj.NOP_rec    = zeros(n1,n2,nFrames);
% Nematic alignment director-x,y dir (eigenvector of nematic order parameter)
orderParamObj.NOPx_rec   = zeros(n1,n2,nFrames);
orderParamObj.NOPy_rec   = zeros(n1,n2,nFrames);
% Averages
orderParamObj.aveC_rec   = zeros(1,nFrames);
orderParamObj.aveP_rec   = zeros(1,nFrames);
orderParamObj.aveN_rec   = zeros(1,nFrames);
% Common parameters
nSqr = n1 * n2;
% loop it
for ii = 1:nFrames
  % calc ops
  [C,POP,POPx,POPy,NOP,NOPx,NOPy] = ...
    OpCPNCalc(n1,n2,Density_rec(:,:,:,ii),...
    phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
  % calculate aveages
  aveC = trapz_periodic( trapz_periodic( real(C), 2), 1 ) ./ nSqr;
  if n3 > 1
    aveP = trapz_periodic( trapz_periodic( real(POP), 2), 1 ) ./ nSqr;
    aveN = trapz_periodic( trapz_periodic( real(NOP), 2), 1 ) ./ nSqr;
  else
    aveP = 0;
    aveN = 0;
  end
  % store it
  orderParamObj.C_rec(:,:,ii) = real(C);
  orderParamObj.POP_rec(:,:,ii)    = real(POP);
  orderParamObj.POPx_rec(:,:,ii) = real(POPx);
  orderParamObj.POPy_rec(:,:,ii) = real(POPy);
  orderParamObj.NOP_rec(:,:,ii)  = real(NOP);
  orderParamObj.NOPx_rec(:,:,ii) = real(NOPx);
  orderParamObj.NOPy_rec(:,:,ii) = real(NOPy);
  orderParamObj.aveC_rec(ii) = aveC;
  orderParamObj.aveP_rec(ii) = aveP;
  orderParamObj.aveN_rec(ii) = aveN;
end

end
