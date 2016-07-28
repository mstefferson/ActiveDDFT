% Make AVIs from Mat
%
% movie files must be the only ones in the folder
function Mat2AviFromList
mkdir('MatMovies')
movefile('Movie*', 'MatMovies')
cd MatMovies
% keyboard
fps     = 2;
% Do it like a newb for now
MovList = ls;
[NumVar, ~] = size(MovList);

for ii = 3:NumVar
    load(MovList(ii,:))
end

movie2avi(MovM_C,'MovM_C.avi','Compression','None','FPS', fps)
movie2avi(MovM_P,'MovM_P.avi','Compression','None','FPS', fps)
movie2avi(MovM_N,'MovM_N.avi','Compression','None','FPS', fps)
movie2avi(MovM_BT,'MovM_BT.avi','Compression','None','FPS', fps)


if 0
    keyboard
    MovList = ls;
    
    [NumVar, ~] = size(MovList);
    
    % Save the name. get rid of blanks
    MovMatStr = strtrim(MovList(6,:));
    % Trim off the .mat
    MovNmStr = MovMatStr(1:end-4);
    load(MovList(6,:));
    aviStr = sprintf('%s.avi',MovNmStr);
    CurrMov = who('-file', MovMatStr);
    movie2avi(MovM_C,aviStr,'Compression','None','FPS', fps)
    for ii = 4:NumVar
        load('MovList(ii,:)')
        movie2avi(mov, 'myPeaks1.avi', 'compression','None', 'fps',10);
        
    end
end

close all
clear 
% keyboard
end


