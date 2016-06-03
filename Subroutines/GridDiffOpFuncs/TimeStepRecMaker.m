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

function [TimeObj] = TimeStepRecMaker(dt,t_tot,t_rec,t_write)

% Stop the user if erroneous parameters are given
if dt > t_rec;
  error('Recorded interval is shorter then timestep fix before proceeding');
end
if t_rec > t_write;
  error('File write interval is shorter than record interval fix before proceeding');
end

% Start building object
TimeObj.dt = dt;

% Fix the recording time to be divisible by the time step
TimeObj.t_rec = floor(t_rec/TimeObj.dt) * TimeObj.dt;


% Fix the write time to be divisible by the record step
TimeObj.t_write = floor(t_write/TimeObj.t_rec) * TimeObj.t_rec;


% Fix the total run time to be divisible by t_write 
TimeObj.t_tot = floor(t_tot/TimeObj.t_write) * TimeObj.t_write;

% Calculate the outputs
% Tot # of time steps
TimeObj.N_time = round(TimeObj.t_tot/TimeObj.dt); 
% Tot # of recorded points. +1 includes 0
TimeObj.N_rec = round(TimeObj.t_tot/TimeObj.t_rec) + 1;
% # of rec points/write chuck
TimeObj.N_recChunk = round(TimeObj.t_write/TimeObj.t_rec );
% # of record chunks / tot time
TimeObj.N_chunks = round(TimeObj.t_tot/TimeObj.t_write); 
% #  time steps / record
TimeObj.N_dtRec = round(TimeObj.t_rec/TimeObj.dt);      
% # of dt/write chuck
TimeObj.N_dtChunk = round(TimeObj.t_write/TimeObj.dt ); 
   
end
