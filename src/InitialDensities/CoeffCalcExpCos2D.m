% CoeffCalcExpLeg2D.m
% Michael Stefferson
%
% Calculates the coefficients for a 2D distribution of the form
% f(\theta) = exp( \sum a_2n cos(2*n*\theta) ) / Z
% Inputs:
% Nc: Number of max coefficients
% x : spatial coordin. vector. x = cos(\theta)
%
% Outputs:
%                   | a_2 0   0 ...  0   |
% CoeffTempMat =    | a_2 a_4 0 ...  0   |
%                   | ...            0
%                   | a_2 a_4 ... a_{2Nc}|

function [Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,phi,bc)
if 1.499 < bc && bc < 1.501
%     keyboard
    error('Initial equilibrium condition does not like bc = 1.5')
end
% Calculate the Legendre coefficients of the kernal
d2nVec   = KernCoeffCalcHardRod2D(Nc);
% Calculate the coefficents of the distribution
CoeffMat      = zeros(Nc,Nc);
% Make the first guess not zero. Iteration gets stuck
CoeffMat(1,1) = rand();
% Condition for while loop. How much coefficient need to be changing by
% before stopping
epsilon = 1e-12;

for i = 1:Nc
    % Make initial guess for the coeff the previous coeff.
    CoeffTempVec_next = CoeffMat(i,:);
    CoeffTempVec      = ones(1,length(CoeffMat(i,:)));
    
    % Iterate the coefficients
    counter = 0;
    while max(abs(CoeffTempVec - CoeffTempVec_next)) > epsilon
        CoeffTempVec = CoeffTempVec_next;
        CosSum = zeros(1, length(phi) );
        % Sum over all the Legendre polynomials
        for n = 1 : i
            CosSum   =  CosSum + CoeffTempVec(n) *  cos(2*n*phi);
        end % End of Legendre poly sums
        
        %         keyboard
        for j = 1:i
            %             keyboard
%             CoeffTempVec_next(j) = pi * bc * d2nVec(j) ...
%                 * ( trapz( phi, cos(2*j*phi) .* exp(CosSum) ) ./ ...
%                 (trapz( phi, exp(CosSum) ) ) );

% If you make phi periodic, need to integrate it as a periodic function
CoeffTempVec_next(j) = pi * bc * d2nVec(j) ...
                * ( trapz_periodic( phi, cos(2*j*phi) .* exp(CosSum) ) ./ ...
                (trapz_periodic( phi, exp(CosSum) ) ) );
        end
        
        counter = counter + 1;
%         keyboard
    end % End coefficient interation
    
%         keyboard
    % Update the matrix
    CoeffMat(i,:)   = CoeffTempVec_next;
    
    % Make a better guess
    if i ~= Nc
        CoeffMat(i+1,:) = CoeffTempVec_next;
    end
    
end % End loop over all coefficients.

Coeff_best = CoeffMat(Nc,:);
% keyboard
end