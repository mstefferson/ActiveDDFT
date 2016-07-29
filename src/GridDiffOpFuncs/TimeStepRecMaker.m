% Function returns parameters used for recording something a distinct time
% intervals and aligns a final run time and recording time
% to fit with a given time step

% Inputs:
% dt: time step size
% t_tot: total run time
% t_rec: record after this much time

% Outputs:
% Nt: Number of time points. Includes zero
% N_rec: Number of recorded time points, including zero
% N_count: Number of time steps to count before recording. 

% Guide: t == timestep
% time    = t*dt
% t_tot   = Nt * dt
% t_rec   = N_count * dt  

function [timeObj] = TimeStepRecMaker(dt,t_tot,t_rec,t_write)

% Fix parameters if erroneous ones are given 
if t_write  > t_tot;
  fprintf('Total time interval is shorter then write interval\n');
  t_write = t_tot;
end

if t_rec > t_write;
  fprintf('File write interval is shorter than record interval. Fix it \n');
  t_rec = t_write;
end

if dt > t_rec;
  fprintf('Recorded interval is shorter then timestep fix before proceeding\n');
  dt = t_rec;
end

% Start building object
timeObj.dt = dt;

% Fix the recording time to be divisible by the time step
timeObj.t_rec = floor(t_rec/timeObj.dt) * timeObj.dt;


% Fix the write time to be divisible by the record step
timeObj.t_write = floor(t_write/timeObj.t_rec) * timeObj.t_rec;


% Fix the total run time to be divisible by t_write 
timeObj.t_tot = floor(t_tot/timeObj.t_write) * timeObj.t_write;

% Calculate the outputs
% Tot # of time steps
timeObj.N_time = round(timeObj.t_tot/timeObj.dt); 
% Tot # of recorded points. +1 includes 0
timeObj.N_rec = round(timeObj.t_tot/timeObj.t_rec) + 1;
% # of rec points/write chuck
timeObj.N_recChunk = round(timeObj.t_write/timeObj.t_rec );
% # of record chunks / tot time
timeObj.N_chunks = round(timeObj.t_tot/timeObj.t_write); 
% #  time steps / record
timeObj.N_dtRec = round(timeObj.t_rec/timeObj.dt);      
% # of dt/write chuck
timeObj.N_dtChunk = round(timeObj.t_write/timeObj.dt ); 
   
end
