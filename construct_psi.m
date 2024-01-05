function phi = construct_psi(input, output, order, delay)

% Constructs regression matrix from measurements
% Input parameters:
%   - input  (N x 1) = exication of the system
%   - output (N x 1) = response of the system
%   - order  (1 x 1) = assumed order of the system
%   - delay  (1 x 1) = number of samples the output is delayed by
% Output:
%   - phi (N-order-delay x 2*order) = regression matrix
% ================================================================

N = length(input);

% Initialize the matrix 'psi' of the correct size
phi = zeros(N - order - delay, 2 * order);

% Insert measurements of the output in 'psi'
for i = 1 : order
    phi(:, i) = - output(order + delay + 1 - i : end - i);
end

% Insert measurements of the input in 'psi'
for i = 1 : order
    phi(:, i + order) = input(order + 1 - i : end - i - delay);
end

end











