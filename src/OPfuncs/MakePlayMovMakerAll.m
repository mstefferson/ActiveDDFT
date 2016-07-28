file = dir('Movie*');
file.name

for i=1:length(file)
   load(file(i).name);
end

movie2avi(M_All,'MovM_All.avi','Compression','None','FPS', 1)
implay('MovM_All.avi',10)
%Change size
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 700 550]);

