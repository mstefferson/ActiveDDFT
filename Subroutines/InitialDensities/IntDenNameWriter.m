% Writes a string for the initial density

function [IntConcStr] =  IntDenNameWriter( IntConc )

  if IntConc == 0
    IntConcStr = 'PlaneWaveIso';
  elseif IntConc == 1
    IntConcStr = 'PlaneWaveEq';
  elseif IntConc == 2
    IntConcStr = 'PlaneWaveNem';
  elseif IntConc == 3
    IntConcStr = 'Loaded';
  elseif IntConc == 4
    IntConcStr = 'SepPWiso';
  elseif IntConc == 5
    IntConcStr = 'SepPWeq';
  elseif IntConc == 6
    IntConcStr = 'Gaussian';
  end
end
