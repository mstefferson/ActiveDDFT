% opGetSnapshot( path2dir, frame )
% Get the ops of a finished run
function opGetSnapshot( path2dir, frame )
% set-up files
opName = dir( [path2dir '/op_*'] );
runName = dir( [path2dir '/run_*'] );
paramName = dir( [path2dir '/params_*'] );
oPFile = matfile([opName.folder '/' opName.name]);
load( [paramName.folder '/' paramName.name] )
runFile = matfile([runName.folder '/' runName.name]);
gridObj = runFile.gridObj;
cScale = particleObj.b;
% store files in ram
oP.C_rec = oPFile.C_rec;
oP.POP_rec = oPFile.POP_rec;
oP.POPx_rec = oPFile.POPx_rec;
oP.POPy_rec = oPFile.POPy_rec;
oP.NOP_rec = oPFile.NOP_rec;
oP.NOPx_rec = oPFile.NOPx_rec;
oP.NOPy_rec = oPFile.NOPy_rec;

% fix frame if too large
maxFrame = length( oPFile.OpTimeRecVec );
frame = min( frame, maxFrame );
oP.C = oP.C_rec(:,:,frame);
oP.POP = oP.POP_rec(:,:,frame);
oP.POPx = oP.POPx_rec(:,:,frame);
oP.POPy = oP.POPy_rec(:,:,frame);
oP.NOP = oP.NOP_rec(:,:,frame);
oP.NOPx = oP.NOPx_rec(:,:,frame);
oP.NOPy = oP.NOPy_rec(:,:,frame);
%
x = gridObj.x1;
y = gridObj.x2;
% plot it
plotOps( oP, x, y, cScale )
