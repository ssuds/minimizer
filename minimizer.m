%% AA 538 HW4 Part B1
% Shreyas Sudhakar
%
% Numerical 1D search for the unconstrained minimum

%% Clearing all previous data in workspace
clear; 

%% (A) Based on initial point alphaL, knowing that function is first locally decreasing at alpha=alphaL, find upper bound alphaU on interval where min f(alpha) resides
fun = @(x)((x-10)^2);
alphaL = 0; % A point where we know the function is locally decreasing

[alphaL,alpha,alphaU,falphaL,falpha,falphaU] = MNBRAK(fun,alphaL); %Function call to bracketing subroutine to find an upper bound, while also returning an alpha value in between and the function evaluated at each of these points
%%
% New Lower Bound alphaL
display(alphaL);
%%
% New Upper Bound alphaU
display(alphaU);
%%
% Point between lower and upper bound, alpha
display(alpha);
%%
% Function evaluated at Lower Bound
display(falphaL);
%%
% Function evaluated at Upper Bound
display(falphaU);
%%
% Function evaluated at alpha between lower and upper bounds
display(falpha);

%% (B) Given an initial interval [alphaL,alphaU] and knowing the function f(alpha) has a minimum in that interval, narrow down the interval to 0.00503 the size of the original interval
tol = 0.00503; % The fraction of the interval that we want to narrow down to
[alphaL,alpha,alphaU] = GOLDEN(alphaL,alphaU,fun,tol); % Narrow down the range through a golden section process defined in function GOLDEN

%%
% New Lower Bound alphaL
display(alphaL);
%%
% New Upper Bound alphaU
display(alphaU);
%%
% Point between lower and upper bound, alpha
display(alpha);

%% (C) Find the minimum of a function f(alpha) over an interval [alphaL,alphaU] by polynomial approximation
[alpha,falpha] = BRENT(alphaL,alphaU,fun,tol);

%%
% Location of minimum
display(alpha);

%%
% Value of function at minimum
display(falpha);