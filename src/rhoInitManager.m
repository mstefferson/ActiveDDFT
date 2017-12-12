function rhoInitObj = rhoInitManager( rhoInit, systemObj, particleObj )
% build structure
rhoInitObj.intCondInput = rhoInit.intCond;
rhoInitObj.type = rhoInit.intCond{1};
% Build equiblibrium density if you can
% Find non-driving steady state
if strcmp( particleObj.type, 'rods') && strcmp( particleObj.interHb, 'mayer' )
  % Number of coefficients
  Nc    = 20;
  % Equilib distribution. Don't let bc = 1.5
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    rhoInitObj.bc = 1.502;
  else
    rhoInitObj.bc = systemObj.bc;
  end
  if systemObj.n3 == 1
    rhoInitObj.feq = [];
  else
  % Calculate coeff
  phi = 2*pi/systemObj.n3 * (0:systemObj.n3-1);
  [Coeff_best,~] = CoeffCalcExpCos2D(Nc,phi,rhoInitObj.bc); 
  % Build equil distribution
  rhoInitObj.feq = DistBuilderExpCos2Dsing(Nc,phi,Coeff_best);
  end
else
  rhoInitObj.feq  = 1 / ( systemObj.l3 ) .* ones( systemObj.n3, 1 );
end

% int cond
if strcmp(rhoInit.intCond{1}, 'iso')
  fprintf('Setting IC to iso\n');
elseif strcmp(rhoInit.intCond{1}, 'nem')
  fprintf('Setting IC to nem\n');
  if length( rhoInit.intCond ) == 1
    fprintf('user error: setting nem shift angle to 0.\n')
    rhoInit.intCond = {'nem',0};
  end
  rhoInitObj.shiftAngle = rhoInit.intCond{2};
elseif strcmp(rhoInit.intCond{1}, 'eq')
  fprintf('Setting IC to eq\n');
elseif strcmp( rhoInit.intCond{1}, 'load')
  fprintf('Setting IC to load\n');
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
  fprintf('Setting IC to delta function in angle\n');
  if length( rhoInit.intCond ) == 1
    fprintf('user error: setting polar shift angle to 0.\n')
    rhoInit.intCond = {'delP',0};
  end
  rhoInitObj.shiftAngle = rhoInit.intCond{2};
elseif strcmp( rhoInit.intCond{1}, 'crys' )
  fprintf('Setting IC to crystal\n');
  if length( rhoInit.intCond ) == 1
    fprintf('user error: setting crystal lattice spacing to 2.25.\n')
    rhoInit.intCond{2} = 2.25;
  end
  % don't make people guess gaussian guess width, but it is an option
  if length( rhoInit.intCond ) == 2
    rhoInit.intCond{3} = rhoInit.intCond{2} ./ systemObj.l1;
  end
  % set guess to delta function
  rhoInitObj.latticeSpc = rhoInit.intCond{2};
  rhoInitObj.sigGuess = rhoInit.intCond{3};
elseif strcmp( rhoInit.intCond{1}, 'gauss' )
  fprintf('Setting IC to gauss\n')
  rhoInitObj = buildGaussIcObj( rhoInit.intCond, systemObj );
  rhoInitObj.intCondInput = rhoInit.intCond;
elseif strcmp( rhoInit.intCond{1}, 'lorenz' )
  fprintf('Setting IC to lorenztian\n');
  rhoInitObj = buildLorenzianIcObj( rhoInit.intCond, systemObj );
  rhoInitObj.intCondInput = rhoInit.intCond;
else
  fprintf('Cannot find initial density. Setting to iso\n' );
  rhoInitObj.type = 'iso';
end

% perturbations
rhoInitObj.perturbIds = [];
if isempty( rhoInit.perturb )
  rhoInitObj.perturbList = {'none'};
  rhoInitObj.perturb = 'none';
  rhoInitObj.numPerturb = 0;
elseif strcmp( rhoInit.perturb{1,1}, 'none' )
  rhoInitObj.perturbList = {'none'};
  rhoInitObj.perturb = 'none';
  rhoInitObj.numPerturb = 0;
else
  rhoInitObj.perturbIds = ['_'];
  % fix if it's not a list of lists.
  if ischar( rhoInit.perturb{1} )
    numPerturb = 1;
    rhoInitTemp = rhoInit.perturb;
    rhoInitObj.perturb = cell( 1, numPerturb );
    rhoInit.perturb = cell( 1, numPerturb );
    rhoInit.perturb{1} = rhoInitTemp;
  end
  numPerturb = length( rhoInit.perturb );
  rhoInitObj.perturb = cell( 1, numPerturb );
  rhoInitObj.perturbList = cell( 1, numPerturb );
  % loop over parameters
  for ii = 1:numPerturb
    perturbTemp = rhoInit.perturb{ii};
    if strcmp( perturbTemp{1}, 'pw' )
      rhoInitObj.perturbIds = [rhoInitObj.perturbIds 'p'];
      rhoInitObj.perturbList{ii} = perturbTemp{1};
      perturbObj = buildPwIcObj( perturbTemp, systemObj );
      rhoInitObj.perturb{ii} = perturbObj;
    elseif strcmp( perturbTemp{1}, 'gauss' )
      rhoInitObj.perturbIds = [rhoInitObj.perturbIds 'g'];
      rhoInitObj.perturbList{ii} = perturbTemp{1};
      perturbObj = buildGaussIcObj( perturbTemp, systemObj );
      rhoInitObj.perturb{ii} = perturbObj;
    elseif strcmp( perturbTemp{1}, 'lorenz')
      rhoInitObj.perturbIds = [rhoInitObj.perturbIds 'l'];
      rhoInitObj.perturbList{ii} = perturbTemp{1};
      perturbObj = buildLorenzianIcObj( perturbTemp, systemObj );
      rhoInitObj.perturb{ii} = perturbObj;
    else
      fprintf('Warning, cannot find requested perturbations\n')
    end
  end % loop
  % get rid of empty
  rhoInitObj.perturb = ...
    rhoInitObj.perturb( ~cellfun('isempty', rhoInitObj.perturb) );
  rhoInitObj.perturbList = ...
    rhoInitObj.perturbList( ~cellfun('isempty', rhoInitObj.perturbList) );
  rhoInitObj.numPerturb = length( rhoInitObj.perturbList );
end % if perturbations
% build str
rhoInitObj.fileStr = [rhoInitObj.type rhoInitObj.perturbIds];
end

%%%% functions
function pwObj = buildPwIcObj( pwIc, systemObj )
pwObj.type = 'pw';
if systemObj.n1 == 1
  pwIc{4} = 0;
end
if systemObj.n2 == 1
  pwIc{5} = 0;
end
if systemObj.n3 == 1
  pwIc{6} = 0;
end
% Don't perturb more more than you are allowed to
if( pwIc{4} >= systemObj.n1 / 2 )
  pwIc{4} = floor(systemObj.n1 / 2) - 1;
end
if( pwIc{5} >= systemObj.n2 / 2 )
  pwIc{5} = floor(systemObj.n2 / 2) - 1;
end
if( pwIc{6} >= systemObj.n3 / 2 )
  pwIc{6} = floor(systemObj.n3 / 2) - 1;
end
pwObj.randFlag = pwIc{2};
pwObj.amp = pwIc{3};
pwObj.numModes1 = pwIc{4};
pwObj.numModes2 = pwIc{5};
pwObj.numModes3 = pwIc{6};
end
function gaussObj = buildGaussIcObj( gaussIc, systemObj )
gaussObj.type = 'gauss';
numSet = length( gaussIc );
if numSet < 2
  ampVal = 1e-3;
  fprintf('user error: no amplitude given. setting to %f.\n', ampVal)
  gaussIc{2} = ampVal;
end
if numSet < 4
  fprintf('user error: setting gaussian pos 1 values [ l1/2, l1/2 ].\n')
  gaussIc{3} = systemObj.l1/2;
  gaussIc{4} = systemObj.l1/2;
end
if numSet < 6
  fprintf('user error: setting gaussian pos 2 values [ l2/2, l2/2 ].\n')
  gaussIc{5} = systemObj.l2/2;
  gaussIc{6} = systemObj.l2/2;
end
%  Make sure variance isn't zero if doing polar
if gaussIc{3} == 0; gaussIc{3} = 1; end
if gaussIc{5} == 0; gaussIc{5} = 1; end
gaussObj.amp = gaussIc{2};
gaussObj.var1 = gaussIc{3};
gaussObj.center1 = gaussIc{4};
gaussObj.var2 = gaussIc{5};
gaussObj.center2 = gaussIc{6};
if systemObj.n3 > 1
  if numSet < 7
    fprintf('user error: setting gaussian pos 2 values [ l3/2, l3/2 ].\n')
    gaussIc{7} = systemObj.l3/2;
    gaussIc{8} = 0;
  end
  if gaussIc{7} == 0; gaussIc{7} = 1; end
  gaussObj.var3 = gaussIc{7};
  gaussObj.center3 = gaussIc{8};
end
end
function lorenzObj = buildLorenzianIcObj( lorenzIc, systemObj )
lorenzObj.type = 'lorenz';
numSet = length( lorenzIc );
if numSet < 2
  ampVal = 1e-3;
  fprintf('user error: no amplitude given. setting to %f.\n', ampVal)
  lorenzIc{2} = ampVal;
end
if numSet < 4
  fprintf('user error: setting lorenzian pos 1 values [ l1/2, l1/2 ].\n')
  lorenzIc{3} = systemObj.l1/2;
  lorenzIc{4} = systemObj.l1/2;
end
if numSet < 6
  fprintf('user error: setting lorenzian pos 2 values [ l2/2, l2/2 ].\n')
  lorenzIc{5} = systemObj.l2/2;
  lorenzIc{6} = systemObj.l2/2;
end
% don't let variance be zero
if lorenzIc{3} == 0; lorenzIc{3} = systemObj.l1/2; end
if lorenzIc{5} == 0; lorenzIc{5} = systemObj.l2/2; end
lorenzObj.amp = lorenzIc{2};
lorenzObj.width1 = lorenzIc{3};
lorenzObj.center1 = lorenzIc{4};
lorenzObj.width2 = lorenzIc{5};
lorenzObj.center2 = lorenzIc{6};
if systemObj.n3 > 1
  if numSet < 7
    fprintf('user error: setting lorenzian pos 2 values [ l1/2, l1/2 ].\n')
    lorenzIc{7} = systemObj.l3/2;
    lorenzIc{8} = 0;
  end
  % don't let variance be zero
  if lorenzIc{3} == 0; lorenzIc{7} = 1; end
  lorenzObj.width3 = lorenzIc{7};
  lorenzObj.center3 = lorenzIc{8};
end
end

