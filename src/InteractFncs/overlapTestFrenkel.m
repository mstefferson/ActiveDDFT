% overlap  = overlapTestFrenkel(deltaX, deltaY, v1X, v1Y, v2X, v2Y, lRod)
%
% description: overlap test from Frenkel's paper (Frenkle et al 'Evidence' 1985)
%
% overlap if
% gi = ( r_{ij} dot v_i ) ^ 2 - ( lRod / 2 ) ^ 2 ( 1 - (v_i dot v_j ) ^ 2 ) < 0
% gj = ( r_{ij} dot v_j ) ^ 2 - ( lRod / 2 ) ^ 2 ( 1 - (v_i dot v_j ) ^ 2 ) < 0
%
% where v_i is unit vector normal to ith particle. r_ij is the distance between particles i,j
%
function [ overlap ] = overlapTestFrenkel(deltaX, deltaY, v1X, v1Y, v2X, v2Y, lRod)
% get sizes of angular input, and allocate for overlap
overlap = zeros( size(v1X) );
epsilon = 1e-12; % give it some wiggle room
% check if it's too far
halfRodDistSq = lRod ^ 2 / 4;
if deltaX^2 + deltaY ^ 2 <= lRod ^ 2 % too far
  % calculate the two terms
  secondTerm = halfRodDistSq .* ( 1 - (  ( v1X .* v2X ) + (v1Y .* v2Y ) ).^2 ); 
  firstTerm1 = ( deltaX .* v1X +  deltaY .* v1Y ) .^ 2;
  firstTerm2 = ( deltaX .* v2X +  deltaY .* v2Y ) .^ 2;
  g1 = firstTerm1 -  secondTerm;
  g2 = firstTerm2 -  secondTerm;
  % Frenkels overlap test
  overlapTest = (g1 <= 0+epsilon) + (g2 <= 0+epsilon);
  overlap( overlapTest == 2) = 1;
end
