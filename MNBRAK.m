function [alphaL,alpha,alphaU,falphaL,falpha,falphaU] = MNBRAK(fun,alphaL,alpha);
% Subroutine MNBRAK is based on FORTRAN code from "Numerical Recipes – The Art of Scientific Computing"
% This will search in the downhill direction(defined by the function as
% evaluated at the initial points alphaL and alpha) and 
% return new points alphaL, alpha, alphaU which bracket the function's 
% minimum. It also evaluates the function at those three points. 

% constants used in the code
GOLD = 1.618034; % The default ratio by which successive intervals are magnified
GLIMIT = 100; % The maximum magnification allowed for a parabolic-fit step
TINY = 1E-20;

alpha = alphaL + TINY; % Compute a point slightly to the right of the point alphaL

falphaL = fun(alphaL);
falpha = fun(alpha);

if falpha > falphaL % Switch the roles of alpha and alphaL if they are not decreasing as expected
    temp = alpha;
    alpha = alphaL;
    alphaL = temp;
    temp = falpha;
    falpha = falphaL;
    falphaL = temp;
end

alphaU = alpha+GOLD*(alpha-alphaL); % First guess for alphaU
falphaU = fun(alphaU); % Evaluate the function at this initial guess

while falpha >= falphaU
    %Compute U by parabolic extrapolation from alphaL,alpha,and alphaU.
    %TINY is used to prevent any possible division by zero
    R=(alpha-alphaL)*(falpha-falphaU);
    Q = (alpha-alphaU)*(falpha-falphaL);
    U = alpha-((alpha-alphaU)*Q-(alpha-alphaL)*R)/(2*((sign(Q-R))*(abs(max(abs(Q-R),TINY))))); 
    ULIM = alpha+GLIMIT*(alphaU-alpha);
    if (alpha-U)*(U-alphaU) > 0 % Parabolic U is between B and C
        FU = fun(U);
        if FU < falphaU % Got a minimum between alpha and alphaU
            alphaL = alpha;
            falphaL = falpha;
            alpha = U;
            falpha = FU;
        elseif FU > falpha % Got a minimum between alphaL and U
            alphaU = U;
            falphaU = FU;
        end
        U = alphaU+GOLD*(alphaU-alpha); % Parabolic fit was no use. Use default magnification.
        FU = fun(U);
    elseif (alphaU-U)*(U-ULIM)>0 % Parabolic fit is between C and its allowed limit
        FU = fun(U);
        if FU < falphaU
            alpha = alphaU;
            alphaU = U;
            U = alphaU+GOLD*(alphaU-alpha);
            falpha = falphaU;
            falphaU = FU;
            FU = fun(U);
        end
    elseif (U-ULIM)*(ULIM-alphaU) >= 0 % Limit parabolic U to maximum allowed value
        U = ULIM;
        FU = fun(U);
    else % Reject parabolic U, use default magnification
        U = alphaU+GOLD*(alphaU-alpha);
        FU = fun(U);
    end
    alphaL = alpha; % Eliminate oldest point and continue
    alpha = alphaU;
    alphaU = U;
    falphaL = falpha;
    falpha = falphaU;
    falphaU = FU;
end
end
    
            
    
    
