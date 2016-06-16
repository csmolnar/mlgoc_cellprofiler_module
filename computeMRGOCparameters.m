function parameters = computeMRGOCparameters(alpha_tilde, lambda_tilde, rstar, rhatstar)
% Calculates the contour, phase field and MRF parameters for a stable
% circle. The circle energy has one minimum in range of positive circle radius (it is rstar).
%
% circle energy is proportional to d^2
% possible change it for any function of d = rstar/rhatstar

parameters = struct('alpha', {}, 'lambda', {}, 'beta', {}, 'D', {}, 'rstar', {}, 'rhatstar', {}, 'd', {}, 'epsilon', {});

alpha_C = alpha_tilde;
lambda_C = lambda_tilde*rstar;
beta_C = beta_tilde(alpha_tilde, lambda_tilde, rstar, rhatstar, 'Marie');

parameters(1).alpha = alpha_C;
parameters(1).lambda = lambda_C;
parameters(1).beta = beta_C;
parameters(1).rstar = rstar;
parameters(1).rhatstar = rhatstar;
parameters(1).d = rstar/rhatstar;
parameters(1).epsilon = parameters(1).d;

alpha_f = (3 / 4) * alpha_C;
beta_f = (1 / 4) * beta_C;
lambda_f = lambda_C; % w=4
D_f = lambda_C; % w=4

parameters(2).alpha = alpha_f;
parameters(2).lambda = lambda_f;
parameters(2).beta = beta_f;
parameters(2).D = D_f;
parameters(2).rstar = rstar;
parameters(2).rhatstar = rhatstar;
parameters(2).d = rstar/rhatstar;
parameters(2).epsilon = parameters(2).d;

alpha_M = (1 / 2) * alpha_C;
beta_M = (1 / 4) * beta_C;
D_M = (pi / 16) * lambda_C;

parameters(3).alpha = alpha_M;
parameters(3).beta = beta_M;
parameters(3).D = D_M;

parameters(3).rstar = rstar;
parameters(3).rhatstar = rhatstar;
parameters(3).d = rstar/rhatstar;
parameters(3).epsilon = parameters(3).d;

lambda_M = -lambda_f/4;
parameters(3).lambda = lambda_M;

end

function out = beta_tilde(alpha_tildes, lambda_tilde, rstar, rhatstar, interactionType)
%betat := proc (rstar, rstarhat, lambdatilde, alphat) options operator, arrow; evalf(2*Pi*(lambdatilde+alphat)*rstarhat/dEQhat(rstarhat, rstar)) end proc
% BETA - compute the beta weight from given parameters [hor].
%   Y  = beta(alphas, lambda_c, r0s, d_min, epsilon, interactionType)
%   alpha_tildes - set of alpha_tilde parameters
%   lambda_tilde - lambda_tilde parameter
%   rstar - set of radius of preferred circle
%   rhatstar - modificator
%   interactionType -
%       'Marie'       - interaction function defined [roc03]
%       'Bessel'      - interaction function with Bessel function
%       'exponential' - exponential interaction function
%
%   [mcsaba]:
%   Multi Radius 'goc' model 

[~, R0s] = meshgrid(alpha_tildes, rstar);

fs = evalF10(rstar, rstar/rhatstar,  rstar/rhatstar, interactionType);

[~, Fs] = meshgrid(alpha_tildes, fs);

out = (lambda_tilde +  alpha_tildes) .* R0s ./ Fs;

end

function vF10=evalF10(r0s, d_mins, epsilons, interactionType)
%EVALF10 - evaluate the F10 function [hor] between 0 and pi.
%   Y = evalF10(r0s, d_mins, epsilon, interactionType)
%   r0s - set of radius
%   d_mins, epsilon - function parameters
%   interactionType -
%       'Marie'       - interaction function defined [roc03]
%       'Bessel'      - interaction function with Bessel function
%       'exponential' - exponential interaction function
%
%   [hor]:
%   Computational draft

%   Copyright
%   $Revision: 1.0 $     $Date: $

vF10 = zeros(size(d_mins));

for i = 1:length(d_mins)
    d_min = d_mins(i);
    epsilon = epsilons(i);
    r0 = r0s(i);
    vF10(i ) = 2 * quad(@integrandF10, 0, pi, 10e-8, [], r0, d_min, epsilon, interactionType);
%     vF10(i ) = 2 * integral(@integrandF10, 0, pi, 10e-8, [], r0, d_min, epsilon, interactionType);
end
end

function I=integrandF10(x, r0, d_min, epsilon, interactionType)
%INTEGRANDF10 - returns the value of the F10 function [hor].
%   Y = integrandF00(x, r0, d_min, epsilon, interactionType)
%   r0 - radius
%   x - variable
%   d_min, epsilon - function parameters
%   interactionType -
%       'Marie'       - interaction function defined [roc03]
%       'Bessel'      - interaction function with Bessel function
%       'exponential' - exponential interaction function
%
%   [hor]:
%   Computational draft

%   Copyright 
%   $Revision: 1.0 $     $Date: $

%% Funtion for F10

s_2=sin(x./2);
s_a_2=abs(s_2);

x0=2*r0*s_a_2;
p=phi(d_min, epsilon, x0, interactionType);
p_p=phi_p(d_min, epsilon, x0, interactionType);

I1=p;
I2=r0*s_a_2.*p_p;

I=r0*cos(x).*(I1 + I2);
end


function f=phi(d_min, epsilon, x, interactionType)
%PHI - interaction functions.
%   Y = phi(d_min, epsilon, x, interactionType) returns the intercation
%   between to points.
%   d_min, epsilon - function parameters
%   x - variable
%   interactionType -
%       'Marie'       - interaction function defined [roc03]
%       'Bessel'      - interaction function with Bessel function
%       'exponential' - exponential interaction function
%       'Peter'       - new try with one sinus wav
%   [roc03]:
%   M. Rochery and I. Jermyn and J. Zerubia: Higher Order Active Contours and
%   their Application to the Detection of Line Networks in Satellite
%   Imagery, 2nd IEEE Intl. Workshop on Variational, Geometric and
%   Level Set Methods (VLSM), 2003 Oct.

%   Copyright
%   $Revision: 1.0 $     $Date: $

f = zeros(size(x));

if (strcmp(interactionType, 'exponential'))
    %% exp(-lambda*z) interaction function
    f=exp (-x / d_min);
    
elseif (strcmp(interactionType, 'Bessel'))
    %% Bessel function
    f=(1/(2 * pi)) * besselk(0, x / d_min);
    
elseif (strcmp(interactionType, 'Peter'))
    %%  simple sinus wave
    for i=1: length(x)
        if (x(i) < 0)
            f(i)=0;
        elseif (x(i)>d_min)
            f(i)=0;
        else
            f(i)=cos((x(i)/d_min) * 2 * pi + pi)  * 0.5 + 0.5;
        end
    end
elseif (strcmp(interactionType, 'Marie'))
    %% Marie's function
    for i=1: length(x)
        if (x(i) < d_min-epsilon)
            f(i)=1;
        elseif (x(i)>d_min+epsilon)
            f(i)=0;
        else
            f(i)=0.5 * (1 - (x(i)-d_min)/epsilon - sin(pi*(x(i)-d_min)/epsilon)/pi);
        end
    end
end

end

function f=phi_p(d_min, epsilon, x, interactionType)
%PHI_P - first derivative of the interaction functions.
%   Y = phi_p(d_min, epsilon, x, interactionType) returns the intercation
%   between to points.
%   d_min, epsilon - function parameters 
%   x - variable
%   interactionType -
%       'Marie'       - interaction function defined [roc03]
%       'Bessel'      - interaction function with Bessel function
%       'exponential' - exponential interaction function
%       'Peter'       - new try with one sinus wav
%
%   [roc03]:
%   M. Rochery and I. Jermyn and J. Zerubia: Higher Order Active Contours and 
%   their Application to the Detection of Line Networks in Satellite
%   Imagery, 2nd IEEE Intl. Workshop on Variational, Geometric and 
%   Level Set Methods (VLSM), 2003 Oct.

%   Copyright 
%   $Revision: 1.0 $     $Date: $

f = zeros(size(x));

if (strcmp(interactionType, 'exponential'))
    %% original
    f=-(1 / d_min).*exp(-x / d_min);
    
elseif (strcmp(interactionType, 'Bessel'))
    %% Bessel function    
    f=(1/(2*pi*d_min))*-besselk(1,  x / d_min);
    
elseif (strcmp(interactionType, 'Peter'))
    %%  simple sinus wave

    for i=1: length(x)
        if (x(i) < 0)
            f(i)=0;
        elseif (x(i)>d_min)
            f(i)=0;
        else
            f(i)=- (1/d_min) * sin ((2*x(i)*pi)/d_min + pi) * pi;
        end
    end     
elseif (strcmp(interactionType, 'Marie'))
    %% Marie's function    
    for i=1: length(x)
        if ((x(i) < d_min - epsilon) || (x(i) > d_min + epsilon))
            f(i)=0;
        else
            f(i)=0.5 * (-1-cos(pi*(x(i)-d_min)/epsilon))/epsilon;
        end
    end               
end

end
