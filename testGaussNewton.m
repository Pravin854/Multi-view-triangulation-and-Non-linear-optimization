%% Script to test the Gauss-Newton optimizer

% We consider a Gaussian function as our test case. The function is of the
% form f(x; a, m, s) = a * exp(-0.5*(x - m) / s^2), where 'a' is the
% amplitude, 'm' is the mean, and 's' is the standard deviation of the
% Gaussian curve.


% Number of observations
numObs = 50;

% Parameters of the Gaussian (amplitude, mean, std dev)
a = 10;
m = 0;
s = 20;

% Draw observations from the Gaussian with the given parameters
x = linspace(-25,25,numObs);
y = generateGaussian(x, a, m, s);

% Initialization
a0 = 10;
m0 = 13;
s0 = 19.12;
% Initial parameter vector
p0 = [a0; m0; s0];


% Other constants for Gauss-Newton (GN)

% Max number of iters
maxIters = 3;
% Tolerance on gradient magnitude, cost, and successive step size (i.e.,
% when to stop)
tolerance = 1e-15;


% Setup the coefficient matrices, Jacobian, gradient, etc.

% Parameter vector
p = p0;
% Jacobian
J = getJacobianOfGaussian(x, p);
% Compute the error in the Gaussian, given an estimate of the parameter
% vector, and the observations.
residual = getResidualGaussian(x, p, y);
% Compute error
err = norm(residual,2);
% Store initial error (for analysis)
initErr = err;
% Gradient
g = J'*residual;
% Coefficient matrix
A = J'*J;
% Store errors over time
errValues = [initErr];
% Number of successful and unsuccessful steps
numSuccess = 0;
numFail = 0;

% Stopping criterion (if the gradient is already too small, stop)
stop = (norm(g,'inf') < tolerance);


% LM iterations

for k = 1:maxIters
    
    % Iterate until a feasible point is found
    if ~stop
        % Solve the normal equations for the LM linear system
        deltap = A \ g;
        % Check if the parameter update is less than tolerance
        if norm(deltap,2) < tolerance * norm(p,2)
            stop = true;
        % If it is not, then compute the updated vector (do not update in
        % the original vector, as we will first determine whether or not it
        % decreases the cost).
        else
            % Updated parameter vector
            p = p - deltap;
            % Updated Jacobian
            J = getJacobianOfGaussian(x, p);
            % Updated residuals
            residual = getResidualGaussian(x, p, y);
            % Error resulting from this update
            newErr = norm(residual,2);
            % Updated A, g
            A = J'*J;
            g = J'*residual;
            
            % For analysis
            errValues = [errValues, newErr];
            
            % Determine if stopping criteria is attained
            stop = (norm(g,'inf') <= tolerance) || (err <= tolerance);
        end
    end
    
end

% Predicted Gaussian (using optimized parameters)
y_hat = generateGaussian(x, p(1), p(2), p(3));

% Analysis of optimization
plot(errValues, 'LineWidth', 2, 'Color', 'g');
xlabel('Iterations');
ylabel('Loss');
title('Error plot - GN');
figure;
scatter(x, y, 'filled', 'g');
hold on;
scatter(x, y_hat, 'filled', 'r');
xlabel('X');
ylabel('Y');
title('Esimated parameters - GN');
legend('Ground-Truth', 'Predicted');
fprintf('Number of total iterations: %d\n', k);
fprintf('Number of successful steps: %d\n', numSuccess);
fprintf('Number of unsuccessful steps: %d\n', numFail);
fprintf('Estimated parameters of the Gaussian: ');
p'
fprintf('True parameters of the Gaussian: ');
[a, m, s]
fprintf('Difference in estimated and true parameters: ');
norm(p-[a;m;s],1)
