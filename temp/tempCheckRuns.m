%%
files = dir('Hr*');

for ii = 1:length(files)
  filename = files(ii).name;
  fprintf('%s\n', filename);
  cd(filename)
  paramsTemp = dir('params_*');
  load(paramsTemp.name)
  disp( denRecObj )
  disp( timeObj.dt )
  plotname = dir('MaxC*.fig');
  open( plotname.name )
  plotname = dir('Amp*.fig');
  open( plotname.name )
  cd ..
end