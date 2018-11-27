%% Script to test the Levenberg-Marquardt optimizer

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
a0 = 15;
m0 = 45;
s0 = 15;
% Initial parameter vector
p0 = [a0; m0; s0];


% Other constants for LM

% Max number of iters
maxIters = 50;
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
% Initialize the damping parameter
dampingCoeff = 10e-3 * max(diag(J'*J));
% dampMomentum = 2;
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
        deltap = (A + dampingCoeff*eye(size(A,1))) \ g;
        % Check if the magnitude (L2 norm) of the parameter update is less than 'tolerance'
        % If yes, stop
        % --------  CODE HERE --------- %
        if(norm(deltap) < tolerance)
        stop = 1;
        
        % If it is not, then compute the updated vector (do not update in
        % the original vector, as we will first determine whether or not it
        % decreases the cost).
        else
            % Updated parameter vector
            pnew = p - deltap;
            % Error resulting from this update
            newErr = norm(getResidualGaussian(x, pnew, y),2);
            
            % Check if the new error is less than the previous error by a
            % margin greater than the tolerance. Only then we are in a 
            % trust region. Then, carry out the parameter update.
            if newErr < err
                errValues = [errValues, newErr];
                if abs(newErr - err) < tolerance
                    stop = true;
                    break;
                else
                    numSuccess = numSuccess + 1;
                    % Update the parameter vector, Jacobian, and the
                    % residuals, etc. (update p, J, residual, err, A, g)
                    
                    p = pnew;
                    J = getJacobianOfGaussian(x, p);
                    residual = getResidualGaussian(x, p, y);
                    err = norm(residual,2);
                    A = J'*J;
                    g = J'*residual;
                    
                    % Determine if stopping criteria is attained
                    stop = (norm(g,'inf') <= tolerance) || (err <= tolerance);
                    % Decrease damping coefficient (since we are in the
                    % trust-region), i.e., divide dampingCoeff by 2.
                    % --- code here --- %
                    
                    dampingCoeff = dampingCoeff / 2;

                end
            % We are not in a trust-region
            else
                % disp('Error increased');
                numFail = numFail + 1;
                % Increase the damping coefficient, i.e., multiply dampingCoeff by 2
                % --- code here --- %
                dampingCoeff = dampingCoeff * 2;
            
            end
        end
    end
    
end

% Predicted Gaussian (using optimized parameters)
y_hat = generateGaussian(x, a, m, s);
% Analysis of optimization
plot(errValues, 'LineWidth', 2, 'Color', 'g');
xlabel('Iterations');
ylabel('Loss');
title('Error plot - LM');
figure;
scatter(x, y, 'filled', 'g');
hold on;
scatter(x, y_hat, 'filled', 'r');
xlabel('X');
ylabel('Y');
title('Esimated parameters - LM');
legend('Ground-Truth', 'Predicted');
fprintf('Number of total iterations: %d', k);
fprintf('Number of successful steps: %d', numSuccess);
fprintf('Number of unsuccessful steps: %d', numFail);
fprintf('Estimated parameters of the Gaussian: ');
p'
fprintf('True parameters of the Gaussian: ');
[a, m, s]
fprintf('Difference in estimated and true parameters: ');
norm(p-[a;m;s],1)
