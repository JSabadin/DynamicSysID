function unwrapped_signal = custom_unwrap(signal, window_size, threshold)
    % Custom unwrap function that considers neighboring samples
    % to account for noise in the signal
    % Input parameters:
    %   - signal (Nx1) = input wrapped signal
    %   - window_size (1x1) = number of neighboring samples to consider
    %   - threshold (1x1) = phase jump threshold for unwrapping
    % Output:
    %   - unwrapped_signal (Nx1) = output unwrapped signal
    % ==========================================================

    % Input validation
    if nargin < 3
        error('Not enough input arguments. Please provide signal, window_size, and threshold.');
    end

    if ~isvector(signal) || ~isnumeric(signal)
        error('Input "signal" must be a numeric vector.');
    end

    if ~isscalar(window_size) || ~isnumeric(window_size) || mod(window_size, 1) ~= 0 || window_size < 1
        error('Input "window_size" must be a positive integer.');
    end

    if ~isscalar(threshold) || ~isnumeric(threshold) || threshold <= 0
        error('Input "threshold" must be a positive number.');
    end

    N = length(signal);
    unwrapped_signal = signal;

    % Calculate the phase differences
    phase_diff = diff(signal);

    % Smooth out the noise in the phase differences using a moving average filter
    filter_coeffs = ones(1, window_size) / window_size;
    filtered_phase_diff = filter(filter_coeffs, 1, phase_diff);

    % Unwrap the signal by correcting phase jumps
    for i = 1 : N - 1
        % Check if the filtered phase difference exceeds the threshold
        if abs(filtered_phase_diff(i)) > threshold
            sign_correction = sign(filtered_phase_diff(i));
            unwrapped_signal(i + 1 : end) = unwrapped_signal(i + 1 : end) - sign_correction * 2 * pi;
        end
    end
end
