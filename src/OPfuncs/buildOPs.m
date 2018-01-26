% build all the order parameters
function opSave2 = buildOPs(saveNameOP, systemObj, gridObj, ...
  denRecObj, timeObj, rhoInit, runSave)
opSave2 = matfile(saveNameOP,'Writable',true);
% Save
if systemObj.n3 > 1
  % Commonly used trig functions
  [~,~,phi3D] = meshgrid(gridObj.x2,gridObj.x1,gridObj.x3);
  cosPhi3d = cos(phi3D);
  sinPhi3d = sin(phi3D);
  cos2Phi3d = cosPhi3d .^ 2;
  sin2Phi3d = sinPhi3d .^ 2;
  cossinPhi3d = cosPhi3d .* sinPhi3d;
else
  cosPhi3d=0;sinPhi3d=0;cos2Phi3d=0;sin2Phi3d=0;cossinPhi3d=0;
end
% Build time rec vector
if denRecObj.DidIBreak == 0
  totRec = length( denRecObj.TimeRecVec );
  opTimeRecVec = denRecObj.TimeRecVec ;
else %Don't incldue the blowed up denesity for movies. They don't like it.
  totRec = length( denRecObj.TimeRecVec ) - 1;
  opTimeRecVec = denRecObj.TimeRecVec(1:end-1) ;
end
opSave2.OpTimeRecVec = opTimeRecVec;
% Set up saving
opSave2.C_rec    = zeros(systemObj.n1, systemObj.n2, 2);
if systemObj.n3 > 1
  % Distribution slice
  % currently hardcode out inset
  sliceRho.plotInset = 0;
  dotx1 = round( systemObj.n1/4 ); sliceRho.dotx1 = dotx1;
  doty1 = round( 3 * systemObj.n2/4 ); sliceRho.doty1 = doty1;
  dotx2 = round( systemObj.n1/2 ); sliceRho.dotx2 = dotx2;
  doty2 = round( systemObj.n2/2 ); sliceRho.doty2 = doty2;
  dotx3 = round( 3*systemObj.n1/4 ); sliceRho.dotx3 = dotx3;
  doty3 = round( systemObj.n2/4 ); sliceRho.doty3 = doty3;
  nFrames = length( opTimeRecVec ); sliceRho.nFrames = nFrames;
  sliceRho.slice1 = reshape( runSave.Den_rec( dotx1, doty1, :, 1:nFrames ), ...
    [ systemObj.n3 nFrames] );
  sliceRho.slice2 = reshape( runSave.Den_rec( dotx2, doty2, :, 1:nFrames ), ...
    [ systemObj.n3 nFrames] );
  sliceRho.slice3 = reshape( runSave.Den_rec( dotx3, doty3, :, 1:nFrames ), ...
    [ systemObj.n3 nFrames] );
  sliceRho.phi = gridObj.x3;
  opSave2.sliceRho = sliceRho;
  opSave2.POP_rec  = zeros(systemObj.n1, systemObj.n2, 2);
  opSave2.POPx_rec = zeros(systemObj.n1, systemObj.n2, 2);
  opSave2.POPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
  opSave2.NOP_rec  = zeros(systemObj.n1, systemObj.n2, 2);
  opSave2.NOPx_rec = zeros(systemObj.n1, systemObj.n2, 2);
  opSave2.NOPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
  opSave2.aveC_rec = zeros(1,2);
  opSave2.aveP_rec = zeros(1,2);
  opSave2.aveN_rec = zeros(1,2);
end
% Break it into chunks
numChunks = timeObj.N_chunks;
sizeChunk = floor( totRec/ numChunks );
if sizeChunk > 0
  numChunks = ceil( totRec/ sizeChunk);
else
  numChunks = 1;
end
for i = 1:numChunks
  if i ~= numChunks
    currInd =  (i-1) * sizeChunk + 1: i * sizeChunk;
  else
    if numChunks == 1
      currInd = 1:totRec;
    else
      currInd = (i-1) * sizeChunk:totRec;
    end
  end
  % Make the records
  [OPObjTemp] = CPNrecMaker(systemObj.n1,systemObj.n2,...
    systemObj.n3, opTimeRecVec(currInd), runSave.Den_rec(:,:,:,currInd) ,...
    gridObj.x3,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d );
  % Save it
  opSave2.C_rec(:,:,currInd) = OPObjTemp.C_rec;
  opSave2.aveC_rec(1,currInd) = OPObjTemp.aveC_rec;
  if systemObj.n3 > 1
    opSave2.POP_rec(:,:,currInd) = OPObjTemp.POP_rec;
    opSave2.POPx_rec(:,:,currInd) = OPObjTemp.POPx_rec;
    opSave2.POPy_rec(:,:,currInd) = OPObjTemp.POPy_rec;
    opSave2.NOP_rec(:,:,currInd) = OPObjTemp.NOP_rec;
    opSave2.NOPx_rec(:,:,currInd) = OPObjTemp.NOPx_rec;
    opSave2.NOPy_rec(:,:,currInd) = OPObjTemp.NOPy_rec;
    opSave2.aveC_rec(1,currInd) = OPObjTemp.aveC_rec;
    opSave2.aveP_rec(1,currInd) = OPObjTemp.aveP_rec;
    opSave2.aveN_rec(1,currInd) = OPObjTemp.aveN_rec;
  end
end % loop over chunks
% fix size
opSave2.C_rec = opSave2.C_rec(:,:,1:totRec);
opSave2.aveC_rec = opSave2.aveC_rec(1,1:totRec);
if systemObj.n3 > 1
  opSave2.POP_rec  = opSave2.POP_rec(:,:,1:totRec);
  opSave2.POPx_rec = opSave2.POPx_rec(:,:,1:totRec);
  opSave2.POPy_rec = opSave2.POPy_rec(:,:,1:totRec);
  opSave2.NOP_rec  = opSave2.NOP_rec(:,:,1:totRec);
  opSave2.NOPx_rec = opSave2.NOPx_rec(:,:,1:totRec);
  opSave2.NOPy_rec = opSave2.NOPy_rec(:,:,1:totRec);
  opSave2.aveP_rec = opSave2.aveP_rec(1,1:totRec);
  opSave2.aveN_rec = opSave2.aveN_rec(1,1:totRec);
end
% Now do it for steady state sol
if systemObj.n3 > 1
  % Calc CPN
  [~,~,~,~,opSave2.NOPeq,~,~] = ...
    OpCPNCalc(1, 1, reshape( rhoInit.feq, [1,1,systemObj.n3] ), ...
    gridObj.x3,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
end
