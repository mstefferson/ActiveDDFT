% Function: mayerFncHr.m
%
% Description: Calculates the Mayer Function for infinity thin hard rods for a hard rod at the
% origin with its orientation angle being zero. Using the COM trapezoid intersection
% method (as seen in HardRod---cpp--project.)

function [mayerFnc,mayerFncFt] = mayerFncHr(n1, n2, n3, l1, l2, Lrod)
% allocate
mayerFnc = zeros(n1,n2,n3);
dx    = l1/n1;
dy    = l2/n2;
dphi  = 2*pi / n3;
epsilon = 0.00001;
TooFar = Lrod * Lrod + epsilon;
LrodHSq = TooFar  / 4;
% Loop over x
for i = 1:n1
  % periodic
  if i <= n1/2  + 1
    xTemp = (i-1) * dx;
  else
    xTemp = ( - n1 + (i-1) ) * dx;
  end
  xTempSq = xTemp * xTemp;
  % loop over y
  for j = 1:n2
    % periodic
    if j <= (n2/2+1)
      yTemp = (j-1) * dy;
    else
      yTemp = ( - n2 + (j-1) ) * dy;
    end
    % total distance
    yTempSq = yTemp * yTemp;
    DistSq = yTempSq + xTempSq;
    for k = 1:n3
      % check if it's close enough
      if( DistSq <= TooFar )
        % calculate phi and handle cases seperately
        phiTemp = (k-1) * dphi;
        % phi = 0,pi
        if( k == 1 || k == n3 / 2 + 1 )
          if (xTempSq <= TooFar) && (yTemp == 0)
            mayerFnc(i,j,k) = -1;
          else
            mayerFnc(i,j,k) = 0;
          end
        end
        % phi = pi/2, 3pi/2
        if( k == n3 / 4 + 1 || k == 3 * n3 / 4 +1 )
          if( (yTempSq <= LrodHSq) && xTempSq <= LrodHSq)
            mayerFnc(i,j,k) = -1;
          else
            mayerFnc(i,j,k) = 0;
          end
        end
        % 0 < phi < pi/2
        if( k < floor( n3 / 4) + 1 )
          yuL = tan( phiTemp ) * ( xTemp + Lrod / 2) + epsilon;
          ylL = tan( phiTemp ) * ( xTemp - Lrod / 2) - epsilon;
          yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp );
          if( (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq) )
            mayerFnc(i,j,k) = - 1;
          else
            mayerFnc(i,j,k) = 0;
          end
        end
        % pi/2  < phi < pi
        if( k > floor( n3 / 4) + 1 && k < floor(n3/2) + 1 )
          yuL = tan( phiTemp ) * ( xTemp - Lrod / 2) + epsilon;
          ylL = tan( phiTemp ) * ( xTemp + Lrod / 2) - epsilon;
          yMaxSq =  LrodHSq *  sin( phiTemp ) * sin( phiTemp );
          if (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq)
            mayerFnc(i,j,k) = - 1;
          else
            mayerFnc(i,j,k) = 0;
          end
        end
        % pi < phi < 3pi/2
        if( k > floor(n3/2) + 1 && k < floor(3*n3 / 4) + 1 )
          yuL = tan( phiTemp ) * ( xTemp + Lrod / 2) + epsilon;
          ylL = tan( phiTemp ) * ( xTemp - Lrod / 2) - epsilon;
          yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp);
          if( (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq) )
            mayerFnc(i,j,k) = - 1;
          else
            mayerFnc(i,j,k) = 0;
          end
        end
        % 3pi/2 < phi < 2pi
        if(  k > floor(3*n3 / 4) + 1 )
          yuL = tan( phiTemp ) * ( xTemp - Lrod / 2) + epsilon;
          ylL = tan( phiTemp ) * ( xTemp + Lrod / 2) - epsilon;
          yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp );
          if( (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq) )
            mayerFnc(i,j,k) = - 1;
          else
            mayerFnc(i,j,k) = 0;
          end
        end
      else
        mayerFnc(i,j,k) = 0;
      end  %%end if dist too far
    end %%end k loop
  end %%end y loop
end %%end x loop
mayerFncFt = fftshift(fftn( mayerFnc ) );
end %% function

