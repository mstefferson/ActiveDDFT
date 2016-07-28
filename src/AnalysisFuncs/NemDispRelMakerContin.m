function [lambda] = NemDispRelMakerContin(GridObj,ParamObj,ampl_record,kx,ky,km,kmOrig,NmOrig,D_pos,D_rot,bc)

FT_orig = reshape(ampl_record(:,1),1,NmOrig);

% epsilon = FT_orig - FT_eq

% Find find equilibrium distribution and corresponding fourier coefficients
addpath C:\Users\MWS\Documents\MATLAB\Research\BG\INtransEq
[~, feql] = EqDistMakerMain2D(bc, 10,NmOrig,0);
feql_FT = fftshift(fft(feql));
%Get rid of 0 imaginary part
feql_FT = real(feql_FT);

% Initialize rho
[rho_eq] = DenEq2Drot(GridObj,ParamObj);
kx_holder = 17;
ky_holder = 17;

FT_rho_eq = fftshift( fftn( rho_eq ) );
FT_rho_eqVec = reshape(  FT_rho_eq(kx_holder,ky_holder,:), 1, 32 );


epsilon = FT_orig - FT_rho_eqVec;

% l+n = m. must to an n to satisfy this. km = -N/2..N/2-1

% keyboard
lambda =  - ( D_pos .* ( kx .^ 2 + ky .^ 2 ) + D_rot .* km .^ 2 );
eps_cond = 1e-3;
for mi = 1:length(km)
    miOrig = find(km(mi) == kmOrig);              % the indice of m in original k vector
    if abs( real( epsilon(miOrig) ) ) > eps_cond  % Need to skip modes with zero amplitude
        
        for li = 1:2:NmOrig-1         
            kn = km(mi) - kmOrig(li);
            %         keyboard
            if kn < min(kmOrig)
                kn = kn + NmOrig;
            elseif kn > max(kmOrig)
                kn = kn - NmOrig;
            end
            %Indice
            ni = find(kn == kmOrig);
            if mod(kmOrig(li),2) == 0
                lambda(mi) = lambda(mi) + epsilon(ni)/epsilon(miOrig) .*...
                    4 * pi * bc * feql_FT(li) / NmOrig .* (...
                    D_rot .* ( kmOrig(li) .* km(mi) ) ./ ( kmOrig(li)^2 - 1 ) );
                %             keyboard
            end
            
            if mod( kn , 2 ) == 0 % n even
                lambda(mi) = lambda(mi) + epsilon(ni)/epsilon(miOrig) .*...
                    4 .* pi .* bc .* feql_FT(li) / NmOrig .* (...
                    ( D_pos .* ( kx .^ 2 + ky .^ 2 ) + D_rot .* km(mi) .* kn ) ...
                    ./ ( kn .^ 2 - 1 ) );
            end
        end % end l loop
    end % end if epsilon big enough
end % end m loop

% keyboard
end % end function
% keyboard

