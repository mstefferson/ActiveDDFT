% Test for the intersept by rotating rod one onto the axis, 
% calculating slope of rotated rods 2, and finding intersection point.

function [ overLap ] = overlapTestRotateIntercept(deltaX, deltaY, theta1, theta2, lRod)
% check if it's too far
maxRodDist = lRod ^ 2 / 4;
% vertical slope test
epsilon = 0.00001;
if deltaX^2 + deltaY ^ 2 > maxRodDist % too far
  fprintf('Rods are too far away\n')
  overLap = 0;
else % close enough for chance 
  fprintf('Rods are close enough for possible intersect\n')
  % rotate relative distance into new frame
  deltaXrot = cos(theta1) .* deltaX + sin(theta1) .* deltaY;
  deltaYrot = -sin(theta1) .* deltaX + cos(theta1) .* deltaY;
  relSlope = tan( theta2 - theta1 );
  fprintf(' deltaXrot = %f \n deltaYrot = %f \n relSlope = %f \n',...
    deltaXrot, deltaYrot, relSlope );
  % check for flat case
  if relSlope == 0 
    fprintf('Parallel slopes\n');
    if deltaY' == 0 && deltaX' < lRod ^2
      fprintf('On same line. Overlapping.\n');
      overLap = 1;
    else
      fprintf('Not on same line. No overlap\n');
      overLap = 0;
    end
  else % find intercept
    % check for vertical slope
    if abs( abs( (theta2 - theta1) ) - pi / 2 ) < epsilon || ...
        abs( abs( (theta2 - theta1) ) - 3 * pi / 2 ) < epsilon
      xInter = deltaXrot;
    else
      xInter = deltaXrot - deltaYrot / relSlope;
    end
    fprintf( 'interept distance from center of mass = %f\n', xInter);
    distRod2 = ( deltaXrot -  xInter ) ^ 2 + deltaYrot ^ 2;
    if xInter ^ 2 < maxRodDist && distRod2 < maxRodDist
      fprintf('Rods overlapping\n')
      overLap = 1;
    else
      fprintf('Rods not overlapping\n')
      overLap = 0;
    end % line intercept
  end % else relSlope == 0
end % too far
