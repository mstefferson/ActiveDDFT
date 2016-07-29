function [ShitIsFucked] = ExpTooLongChecker(ticExptemp,ticExpInt,rhoVec_FT,Nx,Ny,Nm,j_record)
ShitIsFucked = 0;
%Make sure things are taking too long. This is a sign density---> inf
if ticExptemp > 50 * ticExpInt
    fprintf('Forgive me, your grace. Using expv is beginning to take a long ass time\n');
    fprintf('I fear that something has gone wrong\n');
    
    %         keyboard
    rho = real(ifftn(ifftshift(reshape( rhoVec_FT,Nx,Ny,Nm ))));
    if j_record == 1
        Density_rec(:,:,:,2) = rho;
    else
        Density_rec(:,:,:,j_record) = rho;
    end
    ShitIsFucked  = 1;
    %         keyboard
end
