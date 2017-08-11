% [paramMat] = MakeParamMat( ParamObj, rhoInit, flags )
%
% Description: Creates a parameter matrix that is used by RunHardRod

function [paramMat, numRuns] = MakeParamMat( systemObj, particleObj, ...
  runObj, rhoInit, potInds, flags )

% Create Paramater matrix
% paramMat columns: (n1, n2, n3, l1, l2, vD, bc, IC, SM, runID)
% runID vector
runID = runObj.runID + (0:(runObj.numTrial-1) );
% handle all Ns the same and square box
if flags.AllNsSame == 1
  if systemObj.n3 == 1
    Nvec = unique( [systemObj.n1 systemObj.n2] );
  else
    Nvec = unique( [systemObj.n1 systemObj.n2 systemObj.n3] );
  end
  n1 = Nvec;
  n2 = 1;
  n3 = 1;
else
  n1 = systemObj.n1;
  n2 = systemObj.n2;
  n3 = systemObj.n3;
end
if flags.SquareBox == 1
  Lvec = unique( [systemObj.l1 systemObj.l2] );
  l1 = Lvec;
  l2 = 1;
else
  l1 = systemObj.l1;
  l2 = systemObj.l2;
end
% long range interaction parameters
if isempty( particleObj.interLr )
  lrL1id = 1;
  lrL2id = 1;
  lrE1id = 1;
  lrE2id = 1;
else
  lrL1id = 1:length(particleObj.lrLs1);
  lrL2id = 1:length(particleObj.lrLs2);
  lrE1id = 1:length(particleObj.lrEs1);
  lrE2id = 1:length(particleObj.lrEs2);
end
% Create parameter matrix using combvec
paramMat = combvec( n1, n2, n3, l1, l2, particleObj.vD, systemObj.bc, ...
  rhoInit.IntCond, flags.StepMeth, runID, ...
  lrL1id, lrL2id, lrE1id, lrE2id, potInds);
% get number of runs
numRuns = size( paramMat, 2 );% Fix all Ns the same and ls
if flags.AllNsSame
  if systemObj.n3 == 1
    paramMat(2,:) = paramMat(1,:);
    paramMat(3,:) = ones(1, numRuns) ;
  else
    paramMat(2,:) = paramMat(1,:);
    paramMat(3,:) = paramMat(1,:);
  end
end
if flags.SquareBox
  paramMat(5,:) = paramMat(4,:);
end
