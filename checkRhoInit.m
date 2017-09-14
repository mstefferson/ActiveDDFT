function rhoInitObj = checkRhoInit( rhoInit, systemObj )
% build structure
rhoInitObj.intCondInput = rhoInit.intCond;
rhoInitObj.pwInput = rhoInit.perturb;
rhoInitObj.intCond = rhoInit.intCond{1};
% int cond
if strcmp(rhoInit.intCond{1}, 'nem')
  if length( rhoInit.intCond ) ~= 1
    fprintf('user error: setting nem shift angle to 0.\n')
    rhoInit.intCond = {'nem',0};
  end
  rhoInitObj.shiftAngle = rhoInit.intCond{2};
elseif strcmp( rhoInit.intCond{1}, 'load')
  if length( rhoInit.intCond ) == 1
    fprintf('user error: no load name. resetting to iso.\n')
    rhoInit.intCond{1} = 'iso';
    if length( rhoInit.intCond ) == 2
      fprintf('setting default load path to ./src/InitialDensities/SavedRhos/')
      rhoInit.intCond{3} = './src/InitialDensities/SavedRhos/';
    end
    if rhoInit.intCond{3}(end) ~= '/'
      rhoInit.intCond{3} = [rhoInit.intCond{3} '/'];
    end
    rhoInitObj.loadName = rhoInit.intCond{2};
    rhoInitObj.pathName = rhoInit.intCond{3};
  elseif strcmp( rhoInit.intCond{1}, 'delP')
    if length( rhoInit.intCond ) ~= 1
      fprintf('user error: setting polar shift angle to 0.\n')
      rhoInit.intCond = {'nem',0};
    end
    rhoInitObj.shiftAngle = rhoInit.intCond{2};
  elseif strcmp( rhoInit.intCond{1}, 'cry' )
    if length( rhoInit.intCond ) == 1
      fprintf('user error: setting crystal lattice spacing to 2.25.\n')
      rhoInit.intCond{2} = 2.25;
      % don't make people guess gaussian guess width, but it is an option
    elseif length( rhoInit.intCond ) == 1
      rhoInit.intCond{3} = rhoInit.intCond{2} ./ systemObj.l1;
    end
    % set guess to delta function
    rhoInitObj.latticeSpc = rhoInit.intCond{2};
    rhoInitObj.sigGuess = rhoInit.intCond{3};
  elseif strcmp( rhoInit.intCond{1}, 'lorenz')
    numSet = length( rhoInit.intCond );
    if numSet < 3
      fprintf('user error: setting lorenzian pos 1 values [ l1/2, l1/2 ].\n')
      rhoInit.intCond{2} = systemObj.l1/2;
      rhoInit.intCond{3} = systemObj.l1/2;
    end
    if numSet < 5
      fprintf('user error: setting lorenzian pos 2 values [ l2/2, l2/2 ].\n')
      rhoInit.intCond{4} = systemObj.l2/2;
      rhoInit.intCond{5} = systemObj.l2/2;
    end
    rhoInitObj.width1 = rhoInit.intCond{2};
    rhoInitObj.center1 = rhoInit.intCond{3};
    rhoInitObj.width2 = rhoInit.intCond{4};
    rhoInitObj.center2 = rhoInit.intCond{5};
    if systemObj.n3 > 1
      rhoInitObj.n3Flag = 1;
      if numSet < 7
        fprintf('user error: setting lorenzian pos 2 values [ l1/2, l1/2 ].\n')
        rhoInit.intCond{6} = systemObj.l3/2;
        rhoInit.intCond{7} = systemObj.l3/2;
      end
      rhoInitObj.width3 = rhoInit.intCond{4};
      rhoInitObj.center3 = rhoInit.intCond{5};
    else
      rhoInitObj.n3Flag = 0;
    end
  end
end

% perturbations
if strcmp( rhoInit.perturb{1}, 'pw' )
  if systemObj.n1 == 1
    rhoInit.perturb{4} = 0;
  end
  if systemObj.n2 == 1
    rhoInit.perturb{5} = 0;
  end
  if systemObj.n3 == 1
    rhoInit.perturb{6} = 0;
  end
  % Don't perturb more more than you are allowed to
  if( rhoInit.perturb{4} >= systemObj.n1 / 2 )
    rhoInit.perturb{4} = floor(systemObj.n1 / 2) - 1;
  end
  if( rhoInit.perturb{5} >= systemObj.n2 / 2 )
    rhoInit.perturb{5} = floor(systemObj.n2 / 2) - 1;
  end
  if( rhoInit.perturb{6} >= systemObj.n3 / 2 )
    rhoInit.perturb{6} = floor(systemObj.n3 / 2) - 1;
  end
  rhoInitObj.randFlag = rhoInit.perturb{2};
  rhoInitObj.perturbWeight = rhoInit.perturb{3};
  rhoInitObj.numModes1 = rhoInit.perturb{4};
  rhoInitObj.numModes2 = rhoInit.perturb{5};
  rhoInitObj.numModes3 = rhoInit.perturb{6};
end
if strcmp( rhoInit.perturb{1}, 'gauss' )
  if length( rhoInit.perturb ) < 7
    fprintf('user error: setting gaussian perturb to default values.\n')
    rhoInit.perturb{2} = 0;
    rhoInit.perturb{3} = 1e-3;
    rhoInit.perturb{4} = systemObj.l1 / 2;
    rhoInit.perturb{5} = systemObj.l1 / 2;
    rhoInit.perturb{6} = systemObj.l2 / 2;
    rhoInit.perturb{7} = systemObj.l2 / 2;
  end
  %  Make sure variance isn't zero if doing polar
  if rhoInit.perturb{3} == 0; rhoInit.perturb{3} = min(systemObj.l1)/2; end
  if rhoInit.perturb{5} == 0; rhoInit.perturb{5} = min(systemObj.l2)/2; end
  % set obj
  rhoInitObj.homoFlag = rhoInit.perturb{2};
  rhoInitObj.amp = rhoInit.perturb{2};
  rhoInitObj.var1 = rhoInit.perturb{3};
  rhoInitObj.center1 = rhoInit.perturb{4};
  rhoInitObj.var2 = rhoInit.perturb{5};
  rhoInitObj.center2 = rhoInit.perturb{6};
  if systemObj.n3 > 1
    if length( rhoInit.perturb ) < 9
      fprintf('user error: setting gaussian perturb to default values.\n')
      rhoInit.perturb{8} = systemObj.l2 / 2;
      rhoInit.perturb{9} = systemObj.l3 / 2;
    end
    if rhoInit.perturb{7} == 0; rhoInit.perturb{7} = min(systemObj.l3)/2; end
    rhoInitObj.var2 = rhoInit{5};
    rhoInitObj.center2 = rhoInit{6};
  end
end
