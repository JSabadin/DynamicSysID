function theta_SA = stochastic_approximation(input, output, order, delay)

% Construct the regression matrix
psi = construct_psi(input, output, order, delay);

% Set the convergence criteria for the STA estimation loop
tolerance = 1e-15;
theta_error = inf;

% Use LS to calculate the initial theta parameter
% theta_SA = ones(2 * order, 1);
theta_SA = psi \ output(order + delay + 1 : end);

theta_SA = theta_SA *0.8;
% size(theta_SA)
% theta_SA = zeros(2*order,1);

N =  length(input(order + delay + 1:end));

% Iteration counter
k = 1;

% Create a matrix to store the theta_SA values at each iteration
theta_SA_history = [];

while (theta_error > tolerance) & (k< N)
    
     % Store the current value of theta_SA in the history matrix
    theta_SA_history = [theta_SA_history, theta_SA];
    

    % Compute the step size for this iteration
    rho = 1 / k;

    % Compute the update to the parameters using the 
    % current predictor and output variables
    theta_SA_new = theta_SA + rho * psi(k+1, :)' * ...
                   (output(order + delay + k + 1) - psi(k + 1, :) * theta_SA);
   
    % Compute the error between the old and new parameter estimates
    theta_error = norm(theta_SA_new - theta_SA);
    
    % Update the parameter estimate for the next iteration
    theta_SA = theta_SA_new;

    % Increment the iteration counter
    k = k + 1;

   

end

% Plot the parameter changes
figure;
plot(theta_SA_history');
xlabel('Iteration');
ylabel('Parameter value');
title('Parameter Convergence Stohastic aproximation');

% Create legend labels for each parameter
legend_labels = cell(1, size(theta_SA_history, 1));
for i = 1:size(theta_SA_history, 1)
    legend_labels{i} = sprintf('Parameter %d', i);
end

legend(legend_labels, 'Location', 'northeastoutside');

end