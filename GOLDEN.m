function [alphaL,alpha,alphaU] = GOLDEN(alphaL,alphaU,fun,tol);
% Subroutine GOLDEN is based on FORTRAN code from "Numerical Recipes – The Art of Scientific Computing"
% Given a function fun, and given a bracketing pair of points alphaL and
% alphaU such that a function minimum is contained between these two
% points, perform a golden section search for the minimum, isolating it to
% a fractional precision of about tol. The new upper and lower bounds
% alphaU and alphaL are returned, as well as the function minimum alpha.

% Golden Ratios
R = 0.61803399;
C = 1-R;

alpha = (alphaL+alphaU)/2; %Compute a value between the upper and lower bounds

% Evaluate the function at each point
falphaL = fun(alphaL);
falphaU = fun(alphaU);
falpha = fun(alpha);

if falpha > falphaL | falpha > falphaU
    display('Something is wrong with the code in GOLDEN!');
end

% At any given time we will keep track of 4 points: X0, X1, X2, X3
X0 = alphaL;
X3 = alphaU;
if abs(alphaL-alphaU) > abs(alpha-alphaL) % Make X0 to X1 the smaller segment
    X1 = alpha;
    X2 = alpha + C*(alphaU-alpha);
else
    X2 = alpha;
    X1 = alpha - C*(alpha-alphaL);
end

% The initial function evaluations. 
F1 = fun(X1); 
F2 = fun(X2);

while abs(X3-X0) > abs(tol*(alphaU-alphaL))
    if F2 < F1 % One possible outcome, its housekeeping
        X0=X1;
        X1 = X2;
        X2 = R*X1+C*X3;
        F0 = F1;
        F1 = F2;
        F2 = fun(X2); % a new function evaluation
    else % The other possible outcome
        X3 = X2;
        X2 = X1;
        X1 = R*X2+C*X0;
        F3 = F2;
        F2 = F1;
        F1 = fun(X1); % the new function evaluation
    end
end

% We are done, return the best of the two current values as alpha
if F1 < F2 
    alpha = X1;
else
    alpha = X2;
end

% Return the narrowed down upper and lower bounds as alphaU and alphaL
alphaL = X0;
alphaU = X3;
end

    
    