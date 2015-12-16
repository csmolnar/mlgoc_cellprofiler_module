function parameters = computeMRGOCIPMparameters(lambda_tilde, rstar, rhatstar)

% possible change it for any function of d = rstar/rhatstar

parameters = struct('alpha', {}, 'lambda', {}, 'beta', {}, 'D', {}, 'rstar', {}, 'rhatstar', {});

alpha_C = alpha_tilde_inflection(lambda_tilde, rstar, rhatstar, 'Marie');
lambda_C = lambda_tilde*rstar;
beta_C = beta_tilde_inflection(lambda_tilde, rstar, rhatstar, 'Marie');

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

function alpha_tilde = alpha_tilde_inflection(lambda_tilde, rstar, rhatstar, interactionType)

F10 = evalF10(rstar, rstar/rhatstar, rstar/rhatstar, interactionType);
F = evalF(rstar, 0, rstar/rhatstar, rstar/rhatstar, interactionType);

alpha_tilde = (lambda_tilde * rstar .* F') ./ (F10 - rstar .* F');

end

function beta_tilde = beta_tilde_inflection(lambda_tilde, rstar, rhatstar, interactionType)

F10 = evalF10(rstar, rstar/rhatstar, rstar/rhatstar, interactionType);
F = evalF(rstar, 0, rstar/rhatstar, rstar/rhatstar, interactionType);

beta_tilde = lambda_tilde * rstar ./ (F10 - rstar .* F');

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

function vF = evalF(r0s, ks, d_min, epsilon, interactionType)
%EVALF - evaluate the F function [hor] between 0 and pi.
%   Y = evalF(r0s, k, d_mins, epsilon, interactionType)
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

warning off;

vF = zeros(length(r0s), length(ks));

for r=1:length(r0s)
    for k=1:length(ks) 
        r0=r0s(r);
        K=ks(k) / r0;
        % N.B. multiply by 2 as integral is only from 0 to pi not from -pi
        % to pi, of an even function.
        vF(r, k)=2 * quad(@integrandF, 0, pi, [], [], r0, K, d_min, epsilon, interactionType); 
%         vF(r, k)=2 * integral(@integrandF, 0, pi, [], [], r0, K, d_min, epsilon, interactionType);
    end    
end
end

function I = integrandF(x, r0, k, d_min, epsilon, interactionType)
%INTEGRANDF - returns the value of the F function [hor].
%   Y = integrandF(x, r0, k, d_min, epsilon, interactionType)
%   r0 - radius
%   x - variable
%   k - Fourier-frequencies
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

s_a_2=abs(sin(x ./ 2));
s_2=sin(x / 2);
c_2=cos(x / 2);
c=cos(x);
s=sin(x);
x0=2*r0*s_a_2;
a=(c_2 .* c_2)./s_a_2;
p=phi(d_min, epsilon, x0, interactionType);
p_p=phi_p(d_min, epsilon, x0, interactionType);
p_d_p=phi_p_p(d_min,epsilon, x0, interactionType);

%%F20%%%%%%%%%%%%%%%%%%%%
I1 = 0.25 * a .* p_p;
I2 = 0.5 * r0 * s_2 .* s_2 .* p_d_p;
I3 = s_a_2 .* p_p;

F20 = r0 * c .* (I1+I2+I3);
% F20a = F20 without divergent term
% F20a = r0 * c .* (I2 + I3);

%%F21%%%%%%%%%%%%%%%%%%
I1 = p;
I2 = 2 * r0 * s_a_2 .* p_p;
I3 = 0.5 * r0 * a .* p_p;
I4 = r0^2 * s_2.*s_2 .* p_d_p;

F21=c.*(I1+I2-I3+I4);
% F21a = F21 without divergent term
% F21a = c.*(I1 + I2 + I4);

%%F2021%%%%%%%%%%%%%%%%%%%
% Divergent terms from F20 and F21 combined
% F2021 = 0.25 .* r0 .* c .* a .*p_p;

%%F23%%%%%%%%%%%%%%%%%%%
I1 = p;
I2 = r0 * s_a_2 .* p_p;

F23 = s.*(I1+I2);

%%F24%%%%%%%%%%%%%%%%%%%%%
F24=c.*p;

I = 2*F20 + cos(r0*k*x).*F21 - 2*r0*k*sin(r0*k*x).*F23 + r0^2*k^2*cos(r0*k*x).*F24;

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

function f=phi_p_p(d_min, epsilon, x, interactionType)
%PHI_P_P - second derivative of the interaction functions.
%   Y = phi_p_p(d_min, epsilon, x, interactionType) returns the intercation
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
    f=d_min^(-2)*exp(-x / d_min);
    
elseif (strcmp(interactionType, 'Bessel'))
    %% Bessel function    
    f=(d_min^(-2)/(2*pi)) * 0.5 * (besselk(0, x / d_min) + besselk(2, x / d_min) );
    
elseif (strcmp(interactionType, 'Peter'))
    %%  simple sinus wave

    for i=1: length(x)
        if (x(i) < 0)
            f(i) = 0;
        elseif (x(i)>d_min)
            f(i) = 0;
        else
            f(i) = -((2*pi*pi)/(d_min * d_min)) * cos((2*x(i)*pi)/d_min + pi) ;
        end;
    end;        
    
elseif (strcmp(interactionType, 'Marie'))
    %% Marie's function
    for i=1: length(x)
        if (abs(x(i) - d_min) > epsilon)
            f(i)=0;
        else
            f(i)=0.5 * sin(pi*(x(i)-d_min)/epsilon)*pi/(epsilon^2);
        end;
    end;               
end  
end