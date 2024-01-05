function theta_IV = instrumental_variables(input, output, order, delay, sample_time)

% Performs the algorithm of instrumental variables, 
% returns model parameters
% Input parameters:
%   - input  (N x 1) = exication of the system
%   - output (N x 1) = response of the system
%   - order  (1 x 1) = assumed order of the system
%   - delay  (1 x 1) = number of samples the output is delayed by
%   - sample_time (1 x 1) = sampling time of the I/O signals (in seconds)
% Output:
%   - theta_IV (2*order x 1) = model parameters
% ================================================================

% Construct the regression matrix
psi = construct_psi(input, output, order, delay);

% Calculate the LS estimator of the system's parameters
theta_LS = psi \ output(order + delay + 1 : end);

% Initialize the Instrumental Variable (IV) estimator to the LS estimator
theta_IV = theta_LS;

% Set the convergence criteria for the IV estimation loop
tolerance = 1e-6;
theta_error = inf;

% Create a matrix to store the theta_IV values at each iteration
theta_IV_history = [];
i = 1;
while (theta_error > tolerance) && (i < 1000)


     % Store the current value of theta_SA in the history matrix
    theta_IV_history = [theta_IV_history, theta_IV];
    

    % Construct a transfer function G using the current
    % IV estimate of the system's parameters
    a = [1; theta_IV(1 : end/2)]';
    b = theta_IV(end/2+1 : end)';
    G = tf(b, a, sample_time);
    
     % Simulate the output of the system using the input data 
     % and the transfer function G
    time = (0 : length(input) - 1)' * sample_time;
    simulated_output = lsim(G, input, time);
    
    % Construct the new regression matrix W - using simulated output
    W = construct_psi(input, simulated_output, order, delay);
    
    % Update the IV estimate of the system's parameters
    % using the W and Psi matrices
    theta_IV_new = pinv(W' * psi + 0.1 * eye(2 * order)) * W' * output(order + delay + 1 : end);
    
    % Calculate the difference between the old and new estimates
    % of the system's parameters
    theta_error = norm(theta_IV_new - theta_IV);
    
    % Update the IV estimate of the system's parameters
    theta_IV = theta_IV_new;

    i = 1 + i;
end


% Plot the parameter changes
figure;
plot(theta_IV_history', 'o-');
xlabel('Iteration');
ylabel('Parameter value');
title('Parameter Convergence Instrumental Variables');

% Create legend labels for each parameter
legend_labels = cell(1, size(theta_IV_history, 1));
for i = 1:size(theta_IV_history, 1)
    legend_labels{i} = sprintf('Parameter %d', i);
end

legend(legend_labels, 'Location', 'northeastoutside');

end











