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