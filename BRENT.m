function [alpha,falpha] = BRENT(alphaL,alphaU,fun,tol)
% Subroutine BRENT is based on FORTRAN code from "Numerical Recipes – The Art of Scientific Computing"
% Given a function fun, and given bracketing upper and lower bounds alphaU
% and alphaL, this routine isolates the minimum to a fractional precision
% of about tol using Brent's method. The minimum function position is
% returned as alpha, and the minimum function value is returned as falpha.

ITMAX = 100; % Maximum allowed number of iterations
CGOLD = 0.3819660; % Golden Ratio
ZEPS = 1.0E-10; % Small number which protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero

alpha = (alphaL+alphaU)/2; %Compute a value between the upper and lower bounds

% Evaluate the function at each point
falphaL = fun(alphaL);
falphaU = fun(alphaU);
falpha = fun(alpha);

if falpha > falphaL | falpha > falphaU
    display('Something is wrong with the code in GOLDEN!');
end

% A and B must be in ascending order, although the inputs don't have to be
A = min(alphaL,alphaU);
B = max(alphaL,alphaU);

% Initializations
V = alpha;
W = V;
X = V;
E = 0;  % This will be the distance moved on the step before last
FX = fun(X);
FV = FX;
FW = FX;

for iter = 1:ITMAX
    XM = 0.5*(A+B);
    TOL1 = tol*abs(X)+ZEPS;
    TOL2 = 2*TOL1;
    if abs(X-XM) <= TOL2-0.5*(B-A) % Test for done here 
       break 
    end
    
    if abs(E) > TOL1 % Construct a trial parabolic fit
        R = (X-W)*(FX-FV);
        Q = (X-V)*(FX-FW);
        P = (X-V)*Q-(X-W)*R;
        Q = 2*(Q-R);
        if Q > 0
            P = -P;
        end
        Q = abs(Q);
        ETEMP = E;
        E = D;
        if abs(P) >= abs(0.5*Q*ETEMP) | P <= Q*(A-X) | P >= Q*(B-X) % Determine the acceptability of the parabolic fit
            if X >= XM % We arrive here for a golden section step, which we take into the larger of the two segments
                E = A-X;
            else
                E = B-X;
            end
            D = CGOLD*E; % Take the golden section step.
        else % If the parabolic fit is OK
            D = P/Q; % Take the parabolic step
            U = X+D;
            if U-A < TOL2 | B-U < TOL2
                D = sign(XM-X)*abs(TOL1);
            end
        end
    end
    
    if abs(D) >= TOL1 % arrive here with D computed either from parabolic fit, or else from golden section
        U = X+D;
    else
        U = X+(sign(D)*abs(TOL1));
    end
    FU = fun(U); % This is the one function evaluation per iteration, and now we have to decide what to do with our function evaluation
    if FU <= FX
        if U >= X
            A = X;
        else
            B = X;
        end
        V = W;
        FV = FW;
        W = X;
        FW = FX;
        X = U;
        FX = FU;
    else
        if U < X
            A = U
        else
            B = U
        end
        
        if FU <= FW | W == X
            V = W;
            FV = FW;
            W = U;
            FW = FU;
        elseif FU <= FV | V==X | V==W
            V = U;
            FV = FU;
        end
    end
end

% Return the best values
alpha = X;
falpha = FX;

end

        
        
       

