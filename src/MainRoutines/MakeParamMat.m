% [paramMat] = MakeParamMat( ParamObj, rhoInit, flags )
%
% Description: Creates a parameter matrix that is used by RunHardRod

function [paramMat, numRuns] = MakeParamMat( systemObj, particleObj, runObj, rhoInit, flags )
%Find number parameters
numbc = length(systemObj.bc);
numvD = length(particleObj.vD);
numIC = length(rhoInit.IntCond);
numSM = length(flags.StepMeth);
numTr = runObj.num_trial;

% Handls all grid points equal and square box different
if flags.AllNsSame == 1
  Nvec = unique( [systemObj.n1 systemObj.n2 systemObj.n3] );
  numN = length(Nvec);
  numn1 = numN;
  numn2 = 1;
  numn3 = 1;
else
  numn1 = length(systemObj.n1);
  numn2 = length(systemObj.n2);
  numn3 = length(systemObj.n3);
end

if flags.SquareBox == 1
  Lvec = unique( [ParamMaster.l1 ParamMaster.l2] );
  numL = length(Lvec);
  numl1 = numL;
  numl2 = 1;
else
  numl1 = length(systemObj.l1);
  numl2 = length(systemObj.l2);
end

% number of parameters
numn1n2 = numn1*numn2;
numn1n2n3 = numn1n2 * numn3;
numn1n2n3l1 = numn1n2n3 * numl1;
numn1n2n3l1l2 = numn1n2n3l1 * numl2;
numn1n2n3l1l2vD = numn1n2n3l1l2 * numvD;
numParams = numn1n2n3l1l2vD * numbc;
numParamsIC = numParams * numIC;
numParamsICSM = numParamsIC * numSM;
numRuns = numParamsICSM * numTr;

% Create Paramater matrix
% paramMat columns: (n1, n2, n3, l1, l2, vD, bc, IC, SM, runID)
paramMat = zeros( numRuns, 10);

for i = 1:numn1
  for j = 1:numn2
    for k = 1:numn3
      for l = 1:numl1
        for m = 1:numl2
          for n = 1:numvD
            for o = 1:numbc
              for p = 1:numIC
                for q = 1:numSM
                  for r = 1:numTr
                    rowInd =  1 + (i-1) + (j-1) * numn1 + ( k-1) * numn1n2 + ...
                      (l-1) * numn1n2n3 + (m-1) * numn1n2n3l1 + ...
                      (n-1) * numn1n2n3l1l2 + (o-1) * numn1n2n3l1l2vD + ...
                      (p-1) * numParams + (q-1) * numParamsIC +...
                      (r-1) * numParamsICSM;
                    
                    % Special cases
                    if flags.AllNsSame
                      paramMat(rowInd,1:3) = Nvec(i);
                    else
                      paramMat(rowInd,1) = systemObj.n1(i);
                      paramMat(rowInd,2) = systemObj.n2(j);
                      paramMat(rowInd,3) = systemObj.n3(k);
                    end
                    
                    if flags.SquareBox
                      paramMat(rowInd,4:5) = Lvec(l);
                    else
                      paramMat(rowInd,4) = systemObj.l1(l);
                      paramMat(rowInd,5) = systemObj.l2(m);
                    end
                    % Everything else
                    paramMat(rowInd,6) = particleObj.vD(n);
                    paramMat(rowInd,7) = systemObj.bc(o);
                    paramMat(rowInd,8) = rhoInit.IntCond(p);
                    paramMat(rowInd,9) = flags.StepMeth(q);
                    paramMat(rowInd,10) = runObj.runID + (r - 1);
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


