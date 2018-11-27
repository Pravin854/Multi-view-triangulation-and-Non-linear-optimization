function J = getJacobianOfGaussian(x, p)
% GETJACOBIANOFGAUSSIAN  Computes the Jacobian of the Gaussian function.
% Assumes that p is a parameter vector comprising of a, m, s in that order.
% These are the amplitude, mean, and standard deviation of the Gaussian
% respectively.


% Extract individual parameters from parameter vector
a = p(1);
m = p(2);
s = p(3);

% Too cryptic! (Note the transpose in the end, to make it 1 x 3)
% J = sum([exp(-0.5*(x-m).^2/s^2); 
%     ((a*(x-m))/(2*s^2)).*exp(-0.5*(x-m).^2/s^2); 
%     ((2*a*(x-m).^2)/s^2) .* exp(-0.5*(x-m).^2 / s^2)], 2)';
J = [exp(-0.5*(x-m).^2/s^2); 
    ((a*(x-m))/(s^2)).*exp(-0.5*(x-m).^2/s^2); 
    ((a*(x-m).^2)/s^3) .* exp(-0.5*(x-m).^2 / s^2)]';

end
