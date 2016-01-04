% Function: MayerFncDiffBtwPntsCalc.m
%
% Description: Calculates the Mayer Function for infinity thin hard rods for a hard rod at the
% origin with its orientation angle being zero. Using the COM trapezoid intersection
% method (as seen in HardRod---cpp--project.)
 
function [MayerFnc] = MayerFncDiffBtwPntsCalc(Nx, Ny, Nm, Lx, Ly, Lrod)


MayerFnc = zeros(Nx,Ny,Nm);
dx    = Lx/Nx;
dy    = Ly/Ny;
dphi  = 2*pi / Nm;
epsilon = 0.00001;
TooFar = Lrod * Lrod + epsilon;
LrodHSq = TooFar  / 4;

for i = 1:Nx
    if i <= Nx/2  + 1
        xTemp = (i-1) * dx;
    else
        xTemp = ( - Nx + (i-1) ) * dx;
    end
    xTempSq = xTemp * xTemp;
    
    for j = 1:Ny
        
        if j <= (Ny/2+1)
            yTemp = (j-1) * dy;
        else
            yTemp = ( - Ny + (j-1) ) * dy;
        end
        
        yTempSq = yTemp * yTemp;
        DistSq = yTempSq + xTempSq;
        
        for k = 1:Nm
            
            if( DistSq <= TooFar )
                phiTemp = (k-1) * dphi;
                
                % phi = 0,pi
                if( k == 1 || k == Nm / 2 + 1 ) 
                    if (xTempSq <= TooFar) && (yTemp == 0)
                        MayerFnc(i,j,k) = -1;
                    else
                        MayerFnc(i,j,k) = 0;
                    end
                end
                
                % phi = pi/2, 3pi/2
                if( k == Nm / 4 + 1 || k == 3 * Nm / 4 +1 )
                    if( (yTempSq <= LrodHSq) && xTempSq <= LrodHSq)
                        MayerFnc(i,j,k) = -1;
                    else
                        MayerFnc(i,j,k) = 0;
                    end
                end
                
                
                % 0 < phi < pi/2
                if( k < floor( Nm / 4) + 1 )
                    yuL = tan( phiTemp ) * ( xTemp + Lrod / 2) + epsilon;
                    ylL = tan( phiTemp ) * ( xTemp - Lrod / 2) - epsilon;
                    yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp );
                    
                    if( (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq) )
                        MayerFnc(i,j,k) = - 1;
                    else
                        MayerFnc(i,j,k) = 0;
                    end
                end
                
                % pi/2  < phi < pi
                if( k > floor( Nm / 4) + 1 && k < floor(Nm/2) + 1 )
                    yuL = tan( phiTemp ) * ( xTemp - Lrod / 2) + epsilon;
                    ylL = tan( phiTemp ) * ( xTemp + Lrod / 2) - epsilon;
                    yMaxSq =  LrodHSq *  sin( phiTemp ) * sin( phiTemp );
                    
                    if (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq)
                        MayerFnc(i,j,k) = - 1;
                    else
                        MayerFnc(i,j,k) = 0;
                    end
                end
                
                
                % pi < phi < 3pi/2
                if( k > floor(Nm/2) + 1 && k < floor(3*Nm / 4) + 1 )
                    yuL = tan( phiTemp ) * ( xTemp + Lrod / 2) + epsilon;
                    ylL = tan( phiTemp ) * ( xTemp - Lrod / 2) - epsilon;
                    yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp);
                    
                    if( (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq) )
                        MayerFnc(i,j,k) = - 1;
                    else
                        MayerFnc(i,j,k) = 0;
                    end
                end
                
                % 3pi/2 < phi < 2pi
                if(  k > floor(3*Nm / 4) + 1 )
                    yuL = tan( phiTemp ) * ( xTemp - Lrod / 2) + epsilon;
                    ylL = tan( phiTemp ) * ( xTemp + Lrod / 2) - epsilon;
                    yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp );
                    
                    if( (yTemp >= ylL) && (yTemp <= yuL) && (yTempSq <= yMaxSq) )
                        MayerFnc(i,j,k) = - 1;
                    else
                        MayerFnc(i,j,k) = 0;
                    end
                end
                
            else
                MayerFnc(i,j,k) = 0;
            end  %%end if dist too far
        end %%end k loop
    end %%end y loop
end %%end x loop

end %% function

