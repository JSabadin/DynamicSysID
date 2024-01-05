function [theta, P] = recursive_least_squares(u, y, order, delay, lambda)
    % Check if input and output data have the same length
    if length(u) ~= length(y)
        error('Input and output data must have the same length.');
    end

    % Construct the regression matrix
    phi = construct_psi(u, y, order, delay);

    % Initialize variables
    N = size(phi, 1); % number of data points
    n = 2 * order; % number of coefficients to estimate
    theta = zeros(n, 1); % initialize the coefficients vector
    P = 1e3 * eye(n); % initialize the covariance matrix with a large value


    % Create a matrix to store the theta values at each iteration
    theta_history = [];

    % Iterate through the data points
    for k = 1:N
        % Store the current value of theta_SA in the history matrix
         theta_history = [theta_history, theta];

        % Get the current input and output
        phi_k = phi(k, :)';
        y_k = y(k + order + delay);

        % Update covariance matrix
        P = (1 / lambda) * P - (1 / lambda) * P * phi_k * phi_k' * P / (1 + (1 / lambda) * phi_k' * P * phi_k);

        % Update coefficients
        theta = theta + P * phi_k * (y_k - phi_k' * theta);
    end


    % Plot the parameter changes
    figure;
    for i = 1:size(theta_history, 1)
        subplot(size(theta_history, 1), 1, i);
        plot(theta_history(i, :));
        xlabel('Iteration');
        ylabel(sprintf('Param. %d val.', i));
        title(sprintf('Parameter %d Convergence RLS', i));
    end
end
