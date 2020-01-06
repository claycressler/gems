function [x_out] = pick_individuals(x,std,N)

% this function produces a number of traits (=individuals) from a lognormal
% distribution with trait mean x, coefficient of variation in trait cv, and
% number of traits to return N

MU = log(x^2 / sqrt(std^2+x^2)); % mean for lognormal
SIGMA = sqrt(log(std^2/x^2 + 1)); % std for lognormal
x_out = lognrnd(MU,SIGMA,N,1); % pull traits