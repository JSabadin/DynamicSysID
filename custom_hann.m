function [hann_window] = custom_hann(N, Tm)
    % Validate Tm
    if Tm > N
        error('Tm must be less than or equal to N.');
    end

    % Initialize the window
    hann_window = zeros(1, N);

    % Calculate the non-zero interval's start and end indices
    start_index = floor((N - Tm) / 2) + 1;
    end_index = start_index + Tm - 1;

    % Generate the Hann window
    n = 0:(Tm - 1);
    hann_window(start_index:end_index) = 0.5 * (1 - cos(2 * pi * n / (Tm - 1)));
end
