% Function: filelist(strId, DirPath)
% strId = file name identifier
% DirPath = Dir to look in
%
% Description: Returns a cell list of all the file names (strings) in a directory
% 
% Michael Stefferson
% 12-May-2016

function [files] = filelist(strId, DirPath)

if nargin == 0
  strId = '';
  DirPath = pwd;
end

if nargin == 1
  DirPath = pwd;
end

% Get length of file extension for parsing purposes
strLength = length(strId);

% Grab everything in
files = dir( DirPath );
TotItems = length( files );

% Store things we want in a temp then extract what we want;
fileTemp  = cell( TotItems , 1 );
counter = 0;

% Note, converting a sturct field to a cell
for i = 1:TotItems
  if files(i).isdir == 0

    if strLength == 0
      counter  = counter + 1;
      fileTemp{counter} = files(i).name;
    elseif findstr( files(i).name, strId )
      counter  = counter + 1;
      fileTemp{counter} = files(i).name; 
    end
    
  end
end

files  = cell( counter , 1 );
for i = 1:counter
  files(i) = fileTemp(i);
end

end
