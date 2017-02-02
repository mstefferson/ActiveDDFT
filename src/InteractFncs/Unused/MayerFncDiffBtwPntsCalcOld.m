% CURRENTLY NOT WORKING 4-Jan-16
% CURRENTLY NOT USED.  4-Jan-16
% 

% Function: MayerFncDiffBtwPntsCalcOld.m
% 
% Description: Calculates the Mayer Function for infinity thin hard rods using 
% a line intersection method.  Mayer function is calculated for a hard rod at the origin
% with its orientation angle being zero---rod at origin, lying along x-axis.

function [FmDistBtwnPts] = MayerFncDiffBtwPntsCalc(n1, n2, n3, l1, l2, dx, dy, dphi, L_rod)

% Vector of position/ang differences
% Handle even and odd rods differently
if mod(n1,2) == 0 %even    
    xv = [ (0 : dx: l1 / 2) ,  (( -l1/2 + dx ) : dx : -dx)  ]; 
else
    xv = [ (0 : dx: L_box/2 - dx/2) , (-l1 / 2 + dx/2) : dx : -dx ];
end
if mod(n2,2) == 0 %even    
    yv = [ (0 : dy: l2 / 2) ,  (( -l2/2 + dy ) : dy : -dy)  ]; 
else
    yv = [ (0 : dy: l2/2 - dy/2) , (-l2 / 2 + dy/2) : dy : -dy ];
end
if mod(n3,2) == 0 %even    
    phiv = [ (0 : dphi: pi) ,  (( -pi + dphi ) : dphi : -dphi)  ];
else
    phiv = [ (0 : dphi: pi - dphi/2) , (-pi + dphi/2) : dphi : -dphi ];
end

[y3d, x3d, phi3d] = meshgrid(yv, xv, phiv);

%Calculate the max and min x,y coordinates of  every point

Xmax = x3d + L_rod / 2 .* abs(cos(phi3d)) ;
Xmin = x3d - L_rod / 2 .* abs(cos(phi3d)) ;

% keyboard
%Initialze matrices
FmSubSet_Int = zeros(n1,n2,n3);  %All points of different angles that intersect in the right range
% FmSubSet_Int2 = zeros(N,N,N);
FmSubSet_Par = zeros(n1,n2,n3);  %All parallel points that intersect in the right range

epsilon = 0.000000001;  %Prevent things from blowing up

DistBtwn = sqrt( ( y3d).^2  +  (x3d).^2);

%Use the equation of a line to calculate the x-intersept

Numerator =   x3d .* tan( phi3d ) - y3d;
Denominator = tan(phi3d) + epsilon;

XposAtIntsct = Numerator  ./ Denominator;
%             YposAtIntsct = tan(phi3d) .* XposAtIntsct + y3d_temp - x3d_temp.*tan(phi3d);
%Make sure rods aren't too far
tolerance = 0.000001;
Not2Far = DistBtwn <= L_rod + tolerance;

% See if the intersection point is within the length of the
% rod we care about
OvrLpEleSrc = abs(XposAtIntsct) < L_rod/2  + tolerance ;

% See if the intersection point is within the length of the
% test rod
OvrLpEleTestRngAbv = XposAtIntsct >= Xmin;
OvrLpEleTestRngBlw = XposAtIntsct <= Xmax;

%find all the common points (matrix indexing) that are not 0
%Adding all matrices should give a 5 if all criteria is met
SumHolder = OvrLpEleTestRngBlw + OvrLpEleTestRngAbv + OvrLpEleSrc + Not2Far;
%The Mayer fnc for these elements is -1
FmSubSet_Int(SumHolder == 4) = -1;

%Fm currently misses horizontal and verticle rods
AngTol = 0.0000000001;
HorzElements =  abs( sin(phiv) )  <= AngTol ;
VertElements = abs( cos(phiv) ) <= AngTol ;
FmSubSet_Par( abs(xv) < L_rod,1,HorzElements) = -1;
FmSubSet_Par( 1 ,abs(yv) < L_rod,VertElements) = -1;
%                 keyboard
FmHolder = FmSubSet_Par + FmSubSet_Int;
FmDistBtwnPts = FmHolder < 0;
FmDistBtwnPts = -FmDistBtwnPts;
%              keyboard

