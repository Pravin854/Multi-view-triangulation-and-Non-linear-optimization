function residual = getResidualGaussian(x, p, y)
% GETRESIDUALGAUSSIAN  Compute the residual of the Gaussian function, given
% a set of observations and an estimate of parameters.
% x: set of values at which the Gaussian is to be evaluated
% p: the parameter vector for the amplitude, mean, and standard deviation
% of the Gaussian


residual = ((p(1) * exp( (-0.5 * (x - p(2)).^2) / (p(3)^2))) - y)';

end
