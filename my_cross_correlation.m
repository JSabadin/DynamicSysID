% function [cross_correlation] = my_cross_correlation(x, y)
%     % Check if the input signals are row vectors
%     if size(x, 1) > 1
%         x = x';
%     end
%     if size(y, 1) > 1
%         y = y';
%     end
% 
%     % Initialize cross_correlation vector
%     cross_correlation = zeros(1, length(x) + length(y) - 1);
% 
%     % Zero-padding the signals
%     x_padded = [zeros(1, length(y) - 1), x, zeros(1, length(y) - 1)];
%     y_padded = [y, zeros(1, length(x) - 1)];
% 
%     % Compute cross-correlation
%     for i = 1:length(cross_correlation)
%         cross_correlation(i) = sum(x_padded(i:i + length(y) - 1) .* y_padded(1:length(y)));
%     end
% end

function [cross_correlation] = my_cross_correlation(x, y)
    % Ensure the input signals are row vectors
    if size(x, 1) > 1
        x = x';
    end
    if size(y, 1) > 1
        y = y';
    end
    
    % Compute cross-correlation using the convolution function
    cross_correlation = conv(x, fliplr(y));
end
