function rhoInitObj = rhoInitManager( rhoInit, systemObj )
% build structure
rhoInitObj.intCondInput = rhoInit.intCond;
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
  end
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
end

% perturbations
if isempty( rhoInit.perturb )
  rhoInitObj.perturbList = {'none'};
  rhoInitObj.perturb{1}.type = 'none';
elseif strcmp( rhoInit.perturb{1,1}, 'none' )
  rhoInitObj.perturbList = {'none'};
  rhoInitObj.perturb{1}.type = 'none';
else
  numPerturb = size( rhoInit.perturb, 1 );
  rhoInitObj.perturb = cell( numPerturb, 1 );
  rhoInitObj.perturbList = cell( numPerturb, 1 );
  for ii = 1:numPerturb
    perturbTemp = rhoInit.perturb{ii};
    rhoInitObj.perturb{ii}.type = perturbTemp{1};
    rhoInitObj.perturbList{ii} = perturbTemp{1};
    if strcmp( perturbTemp{1}, 'pw' )
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
      if( perturbTemp{4} >= systemObj.n1 / 2 )
        perturbTemp{4} = floor(systemObj.n1 / 2) - 1;
      end
      if( perturbTemp{5} >= systemObj.n2 / 2 )
        perturbTemp{5} = floor(systemObj.n2 / 2) - 1;
      end
      if( perturbTemp{6} >= systemObj.n3 / 2 )
        perturbTemp{6} = floor(systemObj.n3 / 2) - 1;
      end
      rhoInitObj.perturb{ii}.randFlag = perturbTemp{2};
      rhoInitObj.perturb{ii}.perturbWeight = perturbTemp{3};
      rhoInitObj.perturb{ii}.numModes1 = perturbTemp{4};
      rhoInitObj.perturb{ii}.numModes2 = perturbTemp{5};
      rhoInitObj.perturb{ii}.numModes3 = perturbTemp{6};
    end
    if strcmp( perturbTemp{1}, 'gauss' )
      if length( perturbTemp ) < 7
        fprintf('user error: setting gaussian perturb to default values.\n')
        perturbTemp{2} = 0;
        perturbTemp{3} = 1e-3;
        perturbTemp{4} = systemObj.l1 / 2;
        perturbTemp{5} = systemObj.l1 / 2;
        perturbTemp{6} = systemObj.l2 / 2;
        perturbTemp{7} = systemObj.l2 / 2;
      end
      %  Make sure variance isn't zero if doing polar
      if perturbTemp{3} == 0; perturbTemp{3} = min(systemObj.l1)/2; end
      if perturbTemp{5} == 0; perturbTemp{5} = min(systemObj.l2)/2; end
      % set obj
      rhoInitObj.perturb{ii}.homoFlag = perturbTemp{2};
      rhoInitObj.perturb{ii}.amp = perturbTemp{2};
      rhoInitObj.perturb{ii}.var1 = perturbTemp{3};
      rhoInitObj.perturb{ii}.center1 = perturbTemp{4};
      rhoInitObj.perturb{ii}.var2 = perturbTemp{5};
      rhoInitObj.perturb{ii}.center2 = perturbTemp{6};
      if systemObj.n3 > 1
        if length( perturbTemp ) < 9
          fprintf('user error: setting gaussian perturb to default values.\n')
          perturbTemp{8} = systemObj.l2 / 2;
          perturbTemp{9} = systemObj.l3 / 2;
        end
        if perturbTemp{7} == 0; perturbTemp{7} = min(systemObj.l3)/2; end
        rhoInitObj.perturb{ii}.var2 = rhoInit{5};
        rhoInitObj.perturb{ii}.center2 = rhoInit{6};
      end
    end % gauss
    if strcmp( perturbTemp{1}, 'lorenz')
      numSet = length( perturbTemp{:} );
      if numSet < 3
        fprintf('user error: setting lorenzian pos 1 values [ l1/2, l1/2 ].\n')
        perturbTemp{1} = systemObj.l1/2;
        perturbTemp{3} = systemObj.l1/2;
      end
      if numSet < 5
        fprintf('user error: setting lorenzian pos 2 values [ l2/2, l2/2 ].\n')
        perturbTemp{4} = systemObj.l2/2;
        perturbTemp{5} = systemObj.l2/2;
      end
      rhoInitObj.width1 = perturbTemp{2};
      rhoInitObj.center1 = perturbTemp{3};
      rhoInitObj.width2 = perturbTemp{4};
      rhoInitObj.center2 = perturbTemp{5};
      if systemObj.n3 > 1
        rhoInitObj.n3Flag = 1;
        if numSet < 7
          fprintf('user error: setting lorenzian pos 2 values [ l1/2, l1/2 ].\n')
          perturbTemp{6} = systemObj.l3/2;
          perturbTemp{7} = systemObj.l3/2;
        end
        rhoInitObj.perturb{ii}.width3 = perturbTemp{4};
        rhoInitObj.perturb{ii}.center3 = perturbTemp{5};
      else
        rhoInitObj.n3Flag = 0;
      end
    end % lorenztian
  end % loop
end % if perturbations
