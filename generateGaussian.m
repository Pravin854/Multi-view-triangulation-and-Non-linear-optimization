function y = generateGaussian(x, a, m, s)
% GENERATEGAUSSIAN  Generates values for f(x) from a Gaussian function of
% the form y =  f(x) =  a * exp(-0.5 * (x-m)) / s^2), where x is the point
% at which an observation of f(x) is to be made, and a, m, s are the
% amplitude, mean, and standard deviation of the Gaussian respectively.

y = a * exp( (-0.5 * (x - m).^2) / (s^2));

end
