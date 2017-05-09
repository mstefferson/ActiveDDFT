% overlap test from Frenkel's paper
%
% overlap if
% gi = ( r_{ij} dot v_i ) ^ 2 - ( lRod / 2 ) ^ 2 ( 1 - (v_i dot v_j ) ^ 2 ) < 0
% gj = ( r_{ij} dot v_j ) ^ 2 - ( lRod / 2 ) ^ 2 ( 1 - (v_i dot v_j ) ^ 2 ) < 0
%
% where v_i is unit vector normal to ith particle. r_ij is the distance between particles i,j
%
function [ overLap ] = overlapTestFrenkel(deltaX, deltaY, theta1, theta2, lRod)
% check if it's too far
maxRodDist = lRod ^ 2 / 4;
if deltaX^2 + deltaY ^ 2 > maxRodDist % too far
  fprintf('Rods are too far away\n')
  overLap = 0;
else % close enough for chance 
  fprintf('Rods are close enough for possible intersect\n')
  % perpendicular unit vectors. sign doesn't matter so drop it
  v1X = sin(theta1);
  v1Y = cos(theta1);
  v2X = sin(theta2);
  v2Y = cos(theta2);
  % calculate the two terms
  secondTerm = maxRodDist .* ( 1 - (  ( v1X .* v2X ) + (v1Y .* v2Y ) ) ); 
  firstTerm1 = ( deltaX .* v1X ) .^ 2 + ( deltaY .* v1Y ) .^ 2;
  firstTerm2 = ( deltaX .* v2X ) .^ 2 + ( deltaY .* v2Y ) .^ 2;
  g1 = firstTerm1 -  secondTerm;
  g2 = firstTerm2 -  secondTerm;
  fprintf(' firstTerm1 = %f \n firstTerm2 = %f \n secondTerm = %f \n',...
    firstTerm1, firstTerm2, secondTerm );
  fprintf(' g1 = %f \n g2 = %f \n', g1, g2 );
  % Frenkels overlap test
  if (g1 < 0) && (g2 < 0)
    fprintf('Rods overlapping\n')
    overLap = 1;
  else
    fprintf('Rods not overlapping\n')
    overLap = 0;
  end
end
