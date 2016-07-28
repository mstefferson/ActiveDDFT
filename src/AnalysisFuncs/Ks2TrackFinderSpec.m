% k2plotInd: The indices of the k modes you want to plot in relation to the
% original k vector

function [k2plotInd,Nmodes] = Ks2TrackFinderSpec(NumModesMax,ampl_record,TimeRecVec,min_amp)

[N, ~] =  size( ampl_record );

if mod(NumModesMax,2) == 0 
   error('error:evenModes', 'Must use an odd number of modes for now')
end

% The last column might be NaN or Inf. Cut that out
if (isnan(real( ampl_record(1,end) )) ~= 0) || isinf(real( ampl_record(1,end) )) ~= 0
    ampl_record = ampl_record(:,1:end-1);
    TimeRecVec = TimeRecVec(1:end-1);
end
%Only record the amplitudes that are larger than a certain value in the
%first recording. The find function returns which ky elements this
%cooresponds to. Store these in ky2plot

% See if any amplitudes, the beginning, middle, or end are large enough.
k2plottemp1 = find( abs( real( ampl_record(:,1) ) ) > min_amp);
Length1 = length(k2plottemp1);

if mod( length(TimeRecVec), 2 ) == 0
    k2plottemp2 = find( abs( real( ampl_record(:,length(TimeRecVec)/2) ) ) > min_amp);
else
    k2plottemp2 = find( abs( real( ampl_record(:,(length(TimeRecVec) + 1 ) /2) ) ) > min_amp);
end
Length2 = length(k2plottemp2);

k2plottemp3 = find( abs( real( ampl_record(:,end) ) ) > min_amp);
Length3 = length(k2plottemp3);

% Pick the largest one
if Length1 > Length2 && Length1 > Length3
    k2plotInd = k2plottemp1;
elseif Length2 > Length1 && Length2 > Length3
    k2plotInd = k2plottemp2;
elseif Length3 > Length2 && Length3 > Length1
    k2plotInd = k2plottemp3;
else
    k2plotInd = k2plottemp3;
end
% keyboard
%Trim it if it's too long.
% keyboard
if length(k2plotInd) > NumModesMax
    if mod(length(k2plotInd),2) == 0 % Even
    k2plotInd = k2plotInd( ( ( length(k2plotInd) ) / 2 + 1 - (NumModesMax-1) / 2 : ...
                  ( length(k2plotInd) ) / 2 + 1 + (NumModesMax-1)/2 )' );
    else %Odd
    k2plotInd = k2plotInd( ( ( length(k2plotInd)+1 ) / 2 -(NumModesMax-1) / 2 : ...
                      (length(k2plotInd)+1) / 2 + (NumModesMax-1)/2 )' );
    end
end

% Fix it if it's even. Trim off a mode
VecLength = length(k2plotInd);
if mod( VecLength , 2) == 0
k2plotInd = k2plotInd(2:end); 
end
Nmodes = length(k2plotInd);      %Number of modes to study
% keyboard
if Nmodes < 1
   error('You are tracking zero modes. not very interesting') 
end

% keyboard
end