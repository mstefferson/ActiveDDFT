files2convert = dir('*.mat');

for ii = 1:length(files2convert);
  filename = files2convert(ii).name;
  matSave = matfile(filename, 'Writable',true);

  vars = whos('-file',filename );
  if ismember('ParamObj', {vars.name})
    ParamObj = matSave.ParamObj;
    systemObj.Nx = ParamObj.Nx;
    systemObj.Ny = ParamObj.Ny;
    systemObj.Nm = ParamObj.Nm;
    systemObj.Lx = ParamObj.Lx;
    systemObj.Ly = ParamObj.Ly;
    systemObj.Lphi = 2 * pi;
    systemObj.bc = ParamObj.bc;
    
    particleObj.lMaj = ParamObj.L_rod;
    particleObj.vD = ParamObj.vD;

    runObj.runID = ParamObj.runID;
    runObj.num_trial = ParamObj.num_trial;
    runObj.trialID = ParamObj.trialID;

    matSave.systemObj = systemObj;
    matSave.particleObj = particleObj;
    matSave.runObj = runObj;
  end
  if ismember('GridObj', {vars.name})
    GridObj = matSave.GridObj;
    matSave.gridObj = GridObj;
  end
  if ismember('TimeObj', {vars.name})
    TimeObj = matSave.TimeObj;
    matSave.timeObj = TimeObj;
  end
  if ismember('Flags', {vars.name})
    Flags = matSave.Flags;
    matSave.flags = Flags;
  end
  if ismember('RhoInit', {vars.name})
    RhoInit = matSave.RhoInit;
    matSave.rhoInit = RhoInit;
  end
end

