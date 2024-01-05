%% DETERMENING STATIC CHARACTERISTIC 
close all, clear all; clc;
warning('off', 'all');

% Time of simulation, sampling time
simulation_time = 60;
sampling_time = 0.2;
number_of_samples = simulation_time/sampling_time;
% Time array
time = (0:sampling_time:(number_of_samples-1)*sampling_time).';
input_signal.time = time;
input_signal.signals.dimensions = 1;


Y = zeros(size(0 : 0.1 :2.5));
i = 1;
for U = 0 : 0.1 :2.5
    input = ones(1, number_of_samples)*(U);
    input_signal.signals.values = input';
    data = sim("project");
    output = data.output_data.Data;
    Y(i) = mean(output(end-80:end));
    i = i  + 1;
    close all;
    figure
    plot(output)
end
% Rise time = 2.2 - 2.3 seconds
% Lets choose 2  -> 10% of that is sampling time of 0.22
close all;

%%
U = 0 : 0.1 :2.5;
figure
hold on
grid minor
plot(U, Y)
xlabel('Input voltage / V')
ylabel('Output voltage / V')
title('Static characteristic of the pilot plant')
plot(U,Y)


% Define the two points
point1 = [0.8, -1.2];
point2 = [1.2, -2.47];

% Calculate the slope and intercept of the line
slope = (point2(2) - point1(2)) / (point2(1) - point1(1));
intercept = point1(2) - slope * point1(1);

% Create a line that goes through the two points
input_line = linspace(point1(1), point2(1), 100);
output_line = slope * input_line + intercept;

% Plot the line on top of the existing plot
plot(input_line, output_line, 'r', 'LineWidth', 2)


% We go from -1.2 to - 2.47 for the output 
% We go drom 0.8 to 1.2 for the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOING ALL EXPERIMENTS WITH Ts = 0.2 seconds, Tr = 2 seconds
% Creating input singal for ETFE and LSE
% Set parameters
sampling_time_02 = 0.2;
amplitude_min = 0.8;
amplitude_max = 1.2;
prbs_duration = 300; % in seconds
step_duration1 = 50; % in seconds
step_duration2 = 14; % in seconds
simulation_time = prbs_duration + step_duration1 + step_duration2; % in seconds
number_of_samples_02 = simulation_time / sampling_time_02;

% Time array
time = (0:sampling_time_02:(number_of_samples_02-1)*sampling_time_02).';

% Create the step_duration1-second amplitude vector
amplitude_vector = ones(1, (step_duration1 / sampling_time_02)) * amplitude_min;

% Generate the prbs_duration-second PRBS (pseudo-random binary sequence)
prbs_signal = prbs(11)'; % Generate PRBS using the idinput function
prbs_signal = prbs_signal(1:(prbs_duration / sampling_time_02)); % Truncate the PRBS to prbs_duration seconds
prbs_signal = (amplitude_max - amplitude_min) * prbs_signal + amplitude_min; % Scale PRBS to be from 0.8 V to 1.2 V

% Create the 10-second step at the end
step_signal = ones(1, (step_duration2 / sampling_time_02)) * amplitude_max;

% Concatenate the three parts of the signal
input = [amplitude_vector, prbs_signal, step_signal];

% Defining simulation settings
simulation_time = time(end);
sampling_time = sampling_time_02;
input_signal.time = time;
input_signal.signals.dimensions = 1;
input_signal.signals.values = input';

% Checking input signal
figure
plot(time,input)
grid on;
xlabel('t / s', 'FontSize', 14) % Set xlabel font size to 14
ylabel('U / V', 'FontSize', 14) % Set ylabel font size to 14
title('Input signal used for identification', 'FontSize', 16) % Set title font size to 16

% Increase the font size of tick labels on the x and y axes
set(gca, 'FontSize', 12)


data = sim("project");
output = data.output_data.Data;

% Ploting output of simulation
figure
plot(time,output)
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to PRBS')



% Operating point
input_02_oo = input(1);
output_02_oo = mean(output(10/sampling_time_02:45/sampling_time_02));


% Cutting signals and substracting mean
output_02 = output((step_duration1 )/sampling_time_02:end) - output_02_oo;
input_02 = input((step_duration1)/sampling_time_02:end) - input_02_oo;
time_02 = time((step_duration1)/sampling_time_02:end);

% Checking if cutted ok
figure
plot(output_02)
hold on
plot(input_02)

xlabel('samples / n')
ylabel('U / V')
title('Input signal and response of a plant (singals are cutted and mean is substracted)')
legend('Output', 'Input')


%% LS
% Order and delay (in samples) of the system
% The difference between 5 and 6 order is verry low. When summing absoulute
% errors -> order: 1 -> 0.1168
%                  2 -> 0.0533
%                  3 -> 0.0335
%                  4 -> 0.0298
%                  5 -> 0.0294
%                  6 -> 0.0294

order = 4;
delay = 0;

% Define the dependent variable y and regression matrix psi
y = output_02(order + delay + 1 : end);
psi = construct_psi(input_02, output_02, order, delay);

% Solve LS problem
theta = psi \ y;

% Get TF parameters
a_02 = [1; theta(1 : end/2)]';
b_02 = theta(end/2+1 : end)';

% Define the discrete transfer function (also taking into account that the output is delayed)
G_02 = tf(b_02, a_02, sampling_time_02,'OutputDelay', delay)

% Simulate the response
sim_output_02 = lsim(G_02, input_02, time_02) + output_02_oo;
time_training = 0:sampling_time_02:(length(sim_output_02)*sampling_time_02 - sampling_time_02);

figure()
plot(time_training,sim_output_02)
hold on;
plot(time_training,output_02 +  output_02_oo)



xlabel('time / s')
ylabel('U / V')
title('Response of the plant and the model to training signal (singals are cutted and mean is substracted)')
legend('LS - Modeled response', 'Plant response')


% Standard deviations of parameters
error = sim_output_02 - (output_02 +  output_02_oo);
Cov_matrix = 1/(number_of_samples_02 - order + delay)*(error'*error)*pinv(psi'*psi);
std_parameters = sqrt(diag(Cov_matrix));
fprintf("LS std of parameters\n");

for i = 1:order
    fprintf("\nstd(a%d) = %f ", i, std_parameters(i));
end


for i = 1:order
    fprintf("\nstd(b%d) = %f ", i, std_parameters(order + i));
end

G_con_LS = d2c(G_02,'ZOH')


%% % Validation signal
% TODO

% Defining validation singal
% Creating validation singal for LSE
% Set parameters
sampling_time_02 = 0.2;
amplitude_min = 0.8;
amplitude_max = 1.2;
prbs_duration = 30; % in seconds
signal_duration1 = 20; % in seconds
signal_duration2 = 4; % in seconds
signal_duration3 = 20;

simulation_time = prbs_duration + signal_duration1 + signal_duration2 * 4 + signal_duration3; % in seconds
number_of_samples_02 = simulation_time / sampling_time_02;

% Create the signal_duration1-second amplitude vector
amplitude_vector = ones(1, (signal_duration1 / sampling_time_02)) * amplitude_min;

% Generate the prbs_duration-second PRBS (pseudo-random binary sequence)
prbs_signal = prbs(10)'; % Generate PRBS using the idinput function
prbs_signal = prbs_signal(1:(prbs_duration / sampling_time_02)); % Truncate the PRBS to prbs_duration seconds
prbs_signal = (amplitude_max - amplitude_min) * prbs_signal + amplitude_min; % Scale PRBS to be from 0.8 V to 1.2 V


% Create the signal_duration2-second step of 1
step_signal1 = ones(1, (signal_duration2 / sampling_time_02)) * 1;

% Create the signal_duration3-second step of 0.9
step_signal2 = ones(1, (signal_duration2 / sampling_time_02)) * 0.9;

% Create the signal_duration4-second step of 0.1
step_signal3 = ones(1, (signal_duration2 / sampling_time_02)) * 1.1;

% Create the signal_duration5-second step of 0.1
step_signal4 = ones(1, (signal_duration2 / sampling_time_02)) * 0.8;


% Chirp signal
f0 = 0;
f1 = 0.4;
min_amplitude = 0.8; % Minimum amplitude of the chirp signal
max_amplitude = 1.2; % Maximum amplitude of the chirp signal

%time vector for chirp
t = 0:sampling_time_02:(signal_duration3 - sampling_time_02);

chirp_signal = chirp(t, f0, 50, f1, 'linear', -90);

% Normalize the chirp signal to a range of 0 to 1
normalized_chirp = (chirp_signal - min(chirp_signal)) / (max(chirp_signal) - min(chirp_signal));

% Scale the chirp signal to the desired amplitude range
scaled_chirp_signal = min_amplitude + (max_amplitude - min_amplitude) * normalized_chirp;



% Concatenate the three parts of the signal
validation_input = [amplitude_vector,step_signal1,step_signal2,step_signal3,step_signal4,scaled_chirp_signal, prbs_signal];
validation_time = 0:sampling_time_02:(simulation_time - sampling_time_02);

% For the s function!
% figure
% plot(validation_time,validation_input)
% xlabel('time / s')
% ylabel('U / V')
% title('Signal entering s function')
% 

% For the lsim!
% figure
% plot(validation_time,validation_input - amplitude_min)
% xlabel('time / s')
% ylabel('U / V')
% title('Signal entering model')

%% LSE Validation

% Simulate the response
validation_Model_output_02 = lsim(G_02, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
validation_Model_time =  0:sampling_time_02:(length(validation_Model_output_02)*sampling_time_02 - sampling_time_02);

figure
plot(validation_Model_time,validation_Model_output_02)
xlabel('time / s')
ylabel('U / V')
title('Output of model')



simulation_time = validation_time(end);
sampling_time = sampling_time_02;
input_signal.time = validation_time;
input_signal.signals.dimensions = 1;
input_signal.signals.values = validation_input';

data = sim("project");
simulated_validation_output = data.output_data.Data;


simulated_validation_output_oo = mean(simulated_validation_output(signal_duration1/sampling_time_02/2 : signal_duration1/sampling_time_02));

% Ploting output of simulation
hold on;
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to validation signal, performed in the operating point')
legend('LS model response', 'Plant response')


% The input signal in the operating point
figure
plot(validation_Model_time,validation_input(signal_duration1/sampling_time_02:end) - amplitude_min)
xlabel('time / s')
ylabel('U / V')
title('The input signal in the operating point, used for validation')


% Calculating mean of abs error
d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_02));
fprintf("LS: MAE = %.3f\n",d)


%% Instrumental variables
% TODO ###
theta_IV = instrumental_variables(input_02', output_02, 4, 0, sampling_time_02);
a_IV = [1; theta_IV(1 : end/2)]';
b_IV = theta_IV(end/2+1 : end)';
G_IV = tf(b_IV, a_IV, sampling_time_02)
 % Simulate the output of the system using the input data 
 % and the transfer function G
time = (0 : length(input_02) - 1)' * sampling_time_02;
simulated_output = lsim(G_IV, input_02, time);


figure
plot(time, simulated_output)
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to Training signal, performed in the operating point')
legend('IV model response', 'Plant response')
hold on;
plot(time_training,output_02)
legend('IV model response', 'Plant response')



validation_Model_output_IV_02 = lsim(G_IV, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
figure
plot(validation_Model_time, validation_Model_output_IV_02)
% Ploting output of simulation
hold on;
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to validation signal, performed in the operating point')
legend('IV model response', 'Plant response')


% Calculating mean of abs error
d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_IV_02));
fprintf("IV: MAE = %.3f\n",d)
%% RECURSIVE LEAST SQUARES
lambda = 1;
[theta_RLS, P] = recursive_least_squares(input_02, output_02, order, delay, lambda);

a_RLS = [1; theta_RLS(1 : end/2)]';
b_RLS = theta_RLS(end/2+1 : end)';
G_RLS = tf(b_RLS, a_RLS, sampling_time_02)
 % Simulate the output of the system using the input data 
 % and the transfer function G
time = (0 : length(input_02) - 1)' * sampling_time_02;
simulated_output = lsim(G_RLS, input_02, time);


figure
plot(time, simulated_output)
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to Training signal, performed in the operating point')
legend('RLS model response', 'Plant response')
hold on;
plot(time_training,output_02)
legend('RLS model response', 'Plant response')



validation_Model_output_RLS_02 = lsim(G_RLS, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
figure
plot(validation_Model_time, validation_Model_output_RLS_02)
% Ploting output of simulation
hold on;
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to validation signal, performed in the operating point')
legend('RLS model response', 'Plant response')


% Calculating mean of abs error
d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_RLS_02));
fprintf("RLS: MAE = %.3f\n",d)

%% STOHASTIC APROXIMATION
theta_SA = stochastic_approximation(input_02, output_02, order, delay);
a_SA = [1; theta_SA(1 : order)]';
b_SA = theta_SA(order + 1 : 2*order)';
G_SA = tf(b_SA, a_SA, sampling_time_02)

 % Simulate the output of the system using the input data 
simulated_output_SA = lsim(G_SA, input_02, time);

figure
plot(time, simulated_output_SA)
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to Training signal, performed in the operating point')
legend('SA model response', 'Plant response')
hold on;
plot(time_training,output_02)
legend('SA model response', 'Plant response')

validation_Model_output_SA_02 = lsim(G_SA, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
figure
plot(validation_Model_time, validation_Model_output_SA_02)
% Ploting output of simulation
hold on;
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
grid on;
xlabel('t / s')
ylabel('U / V')
title('Response of the plant to validation signal, performed in the operating point')
legend('SA model response', 'Plant response')


% Calculating sum of error
d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_SA_02));
fprintf("SA: MAE = %.3f\n",d)



%% ETFE



% Number of samples
N_02 = length(output_02);

% Frequency resolution
F_02 = 1 / (N_02 * sampling_time_02 );

% Maximum frequency (1/2 of sampling frequency)
f_max_02 = 1 / (2 * sampling_time_02);

% Calculte the FT of input and output signals
input_02_fft = fft(input_02)* sampling_time_02;
output_02_fft = fft(output_02') * sampling_time_02;

% Calculate the ETFE
G_bode_02 = output_02_fft ./ input_02_fft;
G_bode_02_half = G_bode_02(1 : N_02/2);


% Frequencies used for Bode plot
frequencies_02 = 0 : F_02 : f_max_02 - F_02;

% Get frequency response of the model
[amplitude_02, phase_02] = bode(G_02, 2 * pi * frequencies_02);
amplitude_02 = amplitude_02(:);
phase_02 = phase_02(:);



G_bode_02_half_phase = angle(G_bode_02_half);
% G_bode_02_half_phase_unwrapped = unwrap(G_bode_02_half_phase);
G_bode_02_half_phase_unwrapped = custom_unwrap(G_bode_02_half_phase,7,pi);

% Bode plo
figure
subplot(2,1,1)
set(gca, 'XScale', 'log')
hold on
grid on
semilogx(frequencies_02 * 2 * pi, 20 * log10(amplitude_02))
semilogx(frequencies_02 * 2 * pi, 20 * log10(abs(G_bode_02_half)) )
ylabel('|G(j\omega)| / dB')
legend('Model obtained with LS', 'ETFE')
title('Amplitude response')


subplot(2,1,2)
set(gca, 'XScale', 'log')
hold on
grid on
semilogx(frequencies_02 * 2 * pi, phase_02)
semilogx(frequencies_02 * 2 * pi, rad2deg(G_bode_02_half_phase_unwrapped))
xlabel('\omega / rad/s')
ylabel('\theta(G(j\omega)) / deg')
title('Phase response')


%% NOISE MODEL
% TODO

total_time = time(end) + 30;

time_vector = 0:sampling_time:total_time-sampling_time;

simulation_time = total_time;
sampling_time = sampling_time_02;
input_signal.time = time_vector;
input_signal.signals.dimensions = 1;
input_signal.signals.values = 0.8*ones(1, round((simulation_time)/sampling_time_02))';



data = sim("project");
output_used_for_noise = data.output_data.Data;
noise = output_used_for_noise(end-N_02+1:end) - mean(output_used_for_noise(end-N_02+1:end));
noise_time = 0:sampling_time_02:(length(noise)*sampling_time_02 - sampling_time_02);
figure
plot(noise_time,noise)
xlabel('time / s')
ylabel('U / V')
title('Noise signal in the operating point')

% FFT of noise
noise_fft = fft(noise)*sampling_time_02;

% Frequency resolution
F_02 = 1 / (N_02 * sampling_time_02);

% Maximum frequency (1/2 of sampling frequency)
f_max_02 = 1 / (2 * sampling_time_02);

% Frequencies used for Bode plot
frequencies_02 = 0 : F_02 : f_max_02 - F_02;



% Determine the parameters of |N(ω)|
N_0 = 2.2;
N_inf = 0.15;
N_g = (N_0 - N_inf) / sqrt(2);
omega_g = 0.10*2*pi;
f_g = omega_g / (2 * pi);

% Model the noise
abs_noise_model = N_inf + (N_0 - N_inf) ./ ...
                  (sqrt(1 + (frequencies_02 ./ f_g).^2));
abs_noise_model = [abs_noise_model fliplr(abs_noise_model)];
abs_noise_model = abs_noise_model';

% Plot the FT of measured noise and its model
figure
hold on
grid on
% plot(frequencies * 2 * pi, abs(noise_fft(1 : N/2)))
% plot(frequencies * 2 * pi, abs_noise_model(1 : N/2))
plot(frequencies_02, abs(noise_fft(1 : N_02/2)))
plot(frequencies_02, abs_noise_model(1 : N_02/2))
legend('Measured', 'Modeled')
% xlabel('\omega / rad/s')
xlabel('f / Hz')
ylabel('|N(\omega)|')
title('FT of measured noise and its model')

G_noise = abs(noise_fft(1:1 : N_02/2)) ./ abs(input_02_fft((1:1 : N_02/2)))';

% b) Standard deviation of the error due to noise - modeled
G_noise_model = abs_noise_model(1:1 : N_02/2) ./ abs(input_02_fft(1:1 : N_02/2))';


figure
hold on
grid on
plot(frequencies_02, G_noise(1 : N_02/2))
plot(frequencies_02, G_noise_model(1 : N_02/2))
xlabel('f / Hz')
legend('Measured', 'Modeled')
title('Standard deviation of the absolute error of |\DeltaG_n(j\omega)|')


figure
subplot(2, 1, 1)
set(gca, 'XScale', 'log')
hold on
grid on
semilogx(frequencies_02, 20 * log10(abs(G_bode_02_half) - abs(G_noise_model)'))
hold on                              
semilogx(frequencies_02, 20 * log10(abs(G_bode_02_half) + abs(G_noise_model)'))
                             
semilogx(frequencies_02, 20 * log10(abs(G_bode_02_half)));
hold on;
semilogx(frequencies_02, 20 * log10(amplitude_02),':')
hold on;
ylabel('|G(j\omega)| / dB')
legend('Lower bound', 'Upper bound','ETFE','LS','Orthogonal correlation')
title('Bode plot of ETFE and its lower and upper bound')


% Plot phase response 
subplot(2, 1, 2)
set(gca, 'XScale', 'log')
hold on
grid on
% semilogx(frequencies_02, rad2deg(G_bode_02_half_phase_unwrapped) - ...
%                               unwrap(angle(G_noise_model))')
% semilogx(frequencies_02, rad2deg(G_bode_02_half_phase_unwrapped) + ...
%                               unwrap(angle(G_noise_model))')
semilogx(frequencies_02, rad2deg(G_bode_02_half_phase_unwrapped));
hold on;
semilogx(frequencies_02, phase_02,':')
xlabel('f / Hz')
ylabel('arg(G(j\omega)) / deg')



%% Orthogonal correlation
% fg = 0.3

% Define the center and limits
center = 0.125;
lower_limit = center / 10;
upper_limit = center * 10;

% Calculate the logarithms (base 10) of the limits
log_lower_limit = log10(lower_limit);
log_upper_limit = log10(upper_limit);

% Create the logarithmically spaced points
num_points = 11;
harmonic_frequencies = logspace(log_lower_limit, log_upper_limit, num_points);


% Define signal properties
Ts = 0.001; % Time step
T = 100; % Signal duration
num_samples = T / Ts; % Number of samples
t = 0:Ts:(T-Ts); % Time vector

% Initialize matrix to store the signals
input_signals = zeros(num_points, num_samples);
cosine_signals = zeros(num_points, num_samples);
sine_signals = zeros(num_points, num_samples);


% Generate the signals
amplitude = 0.2;
offset = 1;

for i = 1:num_points
    f = harmonic_frequencies(i);
    input_signal_ = amplitude * cos(2 * pi * f * t) + offset;
    cosine_signal = cos(2 * pi * f * t);
    sine_signal = sin(2 * pi * f * t);
    
    input_signals(i, :) = input_signal_;
    cosine_signals(i, :) = cosine_signal;
    sine_signals(i, :) = sine_signal;
end


figure
plot(input_signals(1,:))
hold on;
plot(cosine_signals(1,:))
hold on;
plot(sine_signals(1,:))
legend('input signal', 'cosine used for calculation','sine used for calculation')
title('Plotting signals that will be used for calculation of Orthogonal correlation')
xlabel('t / s')
ylabel('U / V')



simulation_time = t(end);
sampling_time = Ts;
input_signal.time = t';
input_signal.signals.dimensions = 1;

% Initialize output signal
output_signals = zeros(num_points, num_samples);

for i = 1:num_points
    input_signal.signals.values = input_signals(i, :)';
    data = sim("project");
    output_signals(i, :) = data.output_data.Data;
end

%%

% Alghorithm
frequency_response = zeros(length(harmonic_frequencies), 2);
% Loop over the frequencies
for i = 1 : length(harmonic_frequencies)
    U_0 = 0.2;

    % Calculate the period length
    omega = 2 * pi * harmonic_frequencies(i);
    period_length = 1 / harmonic_frequencies(i);


    % Calculate the number of complete periods in the signal
    signal_duration = size(input_signals(i,:),2)*Ts;
    num_periods = floor(signal_duration / period_length);


        
    % Determine the index for cutting the signal
    cut_index = num_periods * period_length / Ts;

    % Trim the output signal and substract the mean
    output_trimmed = output_signals(i,1:cut_index);
   
    % Calculate Re{G(jω)} and Im{G(jω)}
    sine_trimmed = sine_signals(i,1:cut_index);
    cosine_trimmed = cosine_signals(i,1:cut_index);


%     Re_G =  2 * Ts / (U_0 * num_periods * period_length) * sum(output_trimmed .* cosine_trimmed);
%     Im_G = -2 * Ts / (U_0 * num_periods * period_length) * sum(output_trimmed .* sine_trimmed);

    Re_G =  2  / (U_0 ) * mean(output_trimmed .* cosine_trimmed);
    Im_G = -2 / (U_0 ) * mean(output_trimmed .* sine_trimmed);
    
    % Store the frequency response
    frequency_response(i, 1) = sqrt(Re_G ^ 2 + Im_G ^ 2);
    frequency_response(i, 2) = rad2deg(angle(Re_G + 1i * Im_G));
end






figure
subplot(2, 1, 1)
set(gca, 'XScale', 'log')
hold on
grid on
ETFE_plot = semilogx(frequencies_02, 20 * log10(abs(G_bode_02_half)));
scatter_handles = zeros(1, length(harmonic_frequencies));
for i = 1 : length(harmonic_frequencies)
    scatter_handles(i) = scatter(harmonic_frequencies(i), 20 * log10(frequency_response(i, 1)), 'filled', 'k');
end
LS_plot = semilogx(frequencies_02, 20 * log10(abs(amplitude_02)));
ylabel('|G(j\omega)| / dB')
title('Bode plot of G(j\omega)')

legend([ETFE_plot, LS_plot, scatter_handles(1)], 'ETFE', 'LS', 'Orthogonal Correlation');

subplot(2, 1, 2)
set(gca, 'XScale', 'log')
hold on
grid on
semilogx(frequencies_02, rad2deg(G_bode_02_half_phase_unwrapped));
for i = 1 : length(harmonic_frequencies)
    scatter(harmonic_frequencies(i), frequency_response(i, 2), 'filled', 'k')
end
hold on;
semilogx(frequencies_02, phase_02)
ylabel('arg(G(j\omega)) / deg')
xlabel('f / Hz')


%% SPECTRAL ANALYSIS


phi_uu = my_cross_correlation(input_02', input_02');
phi_uy = my_cross_correlation(input_02', output_02');
T_m = length(phi_uu)/8;
w = custom_hann(length(phi_uu), 2*T_m);

% figure
% plot(-250*w)
% hold on;
% plot(phi_uy)
% hold on;
% plot(phi_uu)

% custom_hann = 

Phi_uu = fft(phi_uu .* w);
Phi_uy = fft(phi_uy .* w);


G_spectral_analysis = Phi_uy'./Phi_uu';



% Number of samples
N = length(G_spectral_analysis);

% Frequency resolution
F = 1 / (N * sampling_time_02 );

% Maximum frequency (1/2 of sampling frequency)
f_max = 1 / (2 * sampling_time_02);

% Frequencies used for Bode plot
frequencies = 0 : F : f_max - F;


% Bode plot
figure
subplot(2,1,1)
set(gca, 'XScale', 'log')
hold on
grid on
SA_plot = semilogx(frequencies, 20*log10(abs(G_spectral_analysis(1:end/2))));
LS_plot = semilogx(frequencies_02 , 20 * log10(amplitude_02),'g' );
ETFE = semilogx(frequencies_02 , 20 * log10(abs(G_bode_02_half)),'r:');
for i = 1 : length(harmonic_frequencies)
    scatter_handles(i) = scatter(harmonic_frequencies(i), 20 * log10(frequency_response(i, 1)), 'filled', 'k');
end
ylabel('|G(j\omega)| / dB')
legend([SA_plot, LS_plot,ETFE, scatter_handles(1)], 'Spectral analysis', 'LS','ETFE', 'Orthogonal Correlation');
title('Amplitude response')

G_spectral_analysis_phase = angle(G_spectral_analysis);
G_spectral_analysis_phase_unwraped = custom_unwrap(G_spectral_analysis_phase,3,pi);
% G_spectral_analysis_phase_unwraped(1:6) = G_spectral_analysis_phase_unwraped(1:6) + 2*pi;
% G_bode_02_half_phase_unwrapped(1:2) = G_bode_02_half_phase_unwrapped(1:2) + pi;


subplot(2,1,2)
set(gca, 'XScale', 'log')
hold on
grid on
semilogx(frequencies , rad2deg(G_spectral_analysis_phase_unwraped(1:end/2)))
hold on;;
semilogx(frequencies_02 , (phase_02),'g')
hold on;
% G_bode_02_half_phase_unwrapped(13) = G_bode_02_half_phase_unwrapped(13) - pi;
semilogx(frequencies_02 , rad2deg(G_bode_02_half_phase_unwrapped),'r:')
hold on
for i = 1 : length(harmonic_frequencies)
    scatter(harmonic_frequencies(i), frequency_response(i, 2), 'filled', 'k')
end
xlabel('f / Hz')
ylabel('\theta(G(j\omega)) / deg')
title('Phase response')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNCOMMENT HERE FOR RESULTS WITH Ts = 0.4


% 
% 
% %% DOING ALL EXPERIMENTS WITH Ts = 0.4 seconds, Tr = 2 seconds
% % Creating input singal for ETFE and LSE
% % Set parameters
% sampling_time_02 = 0.4;
% amplitude_min = 0.8;
% amplitude_max = 1.2;
% prbs_duration = 300; % in seconds
% step_duration1 = 50; % in seconds
% step_duration2 = 14; % in seconds
% simulation_time = prbs_duration + step_duration1 + step_duration2; % in seconds
% number_of_samples_02 = simulation_time / sampling_time_02;
% 
% % Time array
% time = (0:sampling_time_02:(number_of_samples_02-1)*sampling_time_02).';
% 
% % Create the step_duration1-second amplitude vector
% amplitude_vector = ones(1, (step_duration1 / sampling_time_02)) * amplitude_min;
% 
% % Generate the prbs_duration-second PRBS (pseudo-random binary sequence)
% prbs_signal = prbs(12)'; % Generate PRBS using the idinput function
% prbs_signal = prbs_signal(1:(prbs_duration / sampling_time_02)); % Truncate the PRBS to prbs_duration seconds
% prbs_signal = (amplitude_max - amplitude_min) * prbs_signal + amplitude_min; % Scale PRBS to be from 0.8 V to 1.2 V
% 
% % Create the 10-second step at the end
% step_signal = ones(1, (step_duration2 / sampling_time_02)) * amplitude_max;
% 
% % Concatenate the three parts of the signal
% input = [amplitude_vector, prbs_signal, step_signal];
% 
% % Defining simulation settings
% simulation_time = time(end);
% sampling_time = sampling_time_02;
% input_signal.time = time;
% input_signal.signals.dimensions = 1;
% input_signal.signals.values = input';
% 
% 
% % Increase the font size of tick labels on the x and y axes
% % set(gca, 'FontSize', 12)
% 
% 
% data = sim("project");
% output = data.output_data.Data;
% 
% 
% 
% 
% % Operating point
% input_02_oo = input(1);
% output_02_oo = mean(output(10/sampling_time_02:45/sampling_time_02));
% 
% 
% % Cutting signals and substracting mean
% output_02 = output((step_duration1 )/sampling_time_02:end) - output_02_oo;
% input_02 = input((step_duration1)/sampling_time_02:end) - input_02_oo;
% time_02 = time((step_duration1)/sampling_time_02:end);
% 
% 
% %% LS
% % Order and delay (in samples) of the system
% % The difference between 5 and 6 order is verry low. When summing absoulute
% % errors -> order: 1 -> 52.5718
% %                  2 -> 36.0761
% %                  3 -> 24.6249
% %                  4 -> 21.1596
% %                  5 -> 19.7132
% %                  6 -> 19.4675
% 
% order =4;
% delay = 0;
% 
% % Define the dependent variable y and regression matrix psi
% y = output_02(order + delay + 1 : end);
% psi = construct_psi(input_02, output_02, order, delay);
% 
% % Solve LS problem
% theta = psi \ y;
% 
% % Get TF parameters
% a_02 = [1; theta(1 : end/2)]';
% b_02 = theta(end/2+1 : end)';
% 
% % Define the discrete transfer function (also taking into account that the output is delayed)
% G_02 = tf(b_02, a_02, sampling_time_02,'OutputDelay', delay)
% 
% % Simulate the response
% sim_output_02 = lsim(G_02, input_02, time_02) + output_02_oo;
% time_training = 0:sampling_time_02:(length(sim_output_02)*sampling_time_02 - sampling_time_02);
% 
% figure()
% plot(time_training,sim_output_02)
% hold on;
% plot(time_training,output_02 +  output_02_oo)
% 
% xlabel('time / s')
% ylabel('U / V')
% title('Response of the plant and the model to training signal (singals are cutted and mean is substracted), Ts = 0.4')
% legend('LS - Modeled response', 'Plant response')
% 
% 
% % Standard deviations of parameters
% error = sim_output_02 - (output_02 +  output_02_oo);
% Cov_matrix = 1/(number_of_samples_02 - order + delay)*(error'*error)*pinv(psi'*psi);
% std_parameters = sqrt(diag(Cov_matrix));
% fprintf("LS std of parameters\n");
% 
% for i = 1:order
%     fprintf("\nstd(a%d) = %f ", i, std_parameters(i));
% end
% 
% 
% for i = 1:order
%     fprintf("\nstd(b%d) = %f ", i, std_parameters(order + i));
% end
% 
% G_con_LS = d2c(G_02,'ZOH')
% 
% 
% %% % Validation signal
% % TODO
% 
% % Defining validation singal
% % Creating validation singal for LSE
% % Set parameters
% sampling_time_02 = 0.4;
% amplitude_min = 0.8;
% amplitude_max = 1.2;
% prbs_duration = 30; % in seconds
% signal_duration1 = 20; % in seconds
% signal_duration2 = 4; % in seconds
% signal_duration3 = 20;
% 
% simulation_time = prbs_duration + signal_duration1 + signal_duration2 * 4 + signal_duration3; % in seconds
% number_of_samples_02 = simulation_time / sampling_time_02;
% 
% % Create the signal_duration1-second amplitude vector
% amplitude_vector = ones(1, (signal_duration1 / sampling_time_02)) * amplitude_min;
% 
% % Generate the prbs_duration-second PRBS (pseudo-random binary sequence)
% prbs_signal = prbs(10)'; % Generate PRBS using the idinput function
% prbs_signal = prbs_signal(1:(prbs_duration / sampling_time_02)); % Truncate the PRBS to prbs_duration seconds
% prbs_signal = (amplitude_max - amplitude_min) * prbs_signal + amplitude_min; % Scale PRBS to be from 0.8 V to 1.2 V
% 
% 
% % Create the signal_duration2-second step of 1
% step_signal1 = ones(1, (signal_duration2 / sampling_time_02)) * 1;
% 
% % Create the signal_duration3-second step of 0.9
% step_signal2 = ones(1, (signal_duration2 / sampling_time_02)) * 0.9;
% 
% % Create the signal_duration4-second step of 0.1
% step_signal3 = ones(1, (signal_duration2 / sampling_time_02)) * 1.1;
% 
% % Create the signal_duration5-second step of 0.1
% step_signal4 = ones(1, (signal_duration2 / sampling_time_02)) * 0.8;
% 
% 
% % Chirp signal
% f0 = 0;
% f1 = 0.4;
% min_amplitude = 0.8; % Minimum amplitude of the chirp signal
% max_amplitude = 1.2; % Maximum amplitude of the chirp signal
% 
% %time vector for chirp
% t = 0:sampling_time_02:(signal_duration3 - sampling_time_02);
% 
% chirp_signal = chirp(t, f0, 50, f1, 'linear', -90);
% 
% % Normalize the chirp signal to a range of 0 to 1
% normalized_chirp = (chirp_signal - min(chirp_signal)) / (max(chirp_signal) - min(chirp_signal));
% 
% % Scale the chirp signal to the desired amplitude range
% scaled_chirp_signal = min_amplitude + (max_amplitude - min_amplitude) * normalized_chirp;
% 
% 
% 
% % Concatenate the three parts of the signal
% validation_input = [amplitude_vector,step_signal1,step_signal2,step_signal3,step_signal4,scaled_chirp_signal, prbs_signal];
% validation_time = 0:sampling_time_02:(simulation_time - sampling_time_02);
% 
% 
% %% LSE Validation
% 
% % Simulate the response
% validation_Model_output_02 = lsim(G_02, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
% validation_Model_time =  0:sampling_time_02:(length(validation_Model_output_02)*sampling_time_02 - sampling_time_02);
% 
% figure
% plot(validation_Model_time,validation_Model_output_02)
% xlabel('time / s')
% ylabel('U / V')
% title('Output of model')
% 
% 
% 
% simulation_time = validation_time(end);
% sampling_time = sampling_time_02;
% input_signal.time = validation_time;
% input_signal.signals.dimensions = 1;
% input_signal.signals.values = validation_input';
% 
% data = sim("project");
% simulated_validation_output = data.output_data.Data;
% 
% 
% simulated_validation_output_oo = mean(simulated_validation_output(signal_duration1/sampling_time_02/2 : signal_duration1/sampling_time_02));
% 
% % Ploting output of simulation
% hold on;
% % plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to validation signal, performed in the operating point, Ts = 0.4 s')
% legend('LS model response', 'Plant response')
% 
% 
% % Calculating mean of abs error
% d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_02));
% fprintf("LS: MAE = %.3f; Ts = 0.4 s \n",d)
% 
% 
% %% Instrumental variables
% % TODO ###
% theta_IV = instrumental_variables(input_02', output_02, 4, 0, sampling_time_02);
% a_IV = [1; theta_IV(1 : end/2)]';
% b_IV = theta_IV(end/2+1 : end)';
% G_IV = tf(b_IV, a_IV, sampling_time_02)
%  % Simulate the output of the system using the input data 
%  % and the transfer function G
% time = (0 : length(input_02) - 1)' * sampling_time_02;
% simulated_output = lsim(G_IV, input_02, time);
% 
% 
% figure
% plot(time, simulated_output)
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to Training signal, performed in the operating point, Ts = 0.4')
% legend('IV model response', 'Plant response')
% hold on;
% plot(time_training,output_02)
% legend('IV model response', 'Plant response')
% 
% 
% 
% validation_Model_output_IV_02 = lsim(G_IV, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
% figure
% plot(validation_Model_time, validation_Model_output_IV_02)
% % Ploting output of simulation
% hold on;
% % plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to validation signal, performed in the operating point, Ts = 0.4')
% legend('IV model response', 'Plant response')
% 
% 
% % Calculating mean of abs error
% d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_IV_02));
% fprintf("IV: MAE = %.3f; Ts = 0.4 s \n",d)
% %% RECURSIVE LEAST SQUARES
% lambda = 1;
% [theta_RLS, P] = recursive_least_squares(input_02, output_02, order, delay, lambda);
% 
% a_RLS = [1; theta_RLS(1 : end/2)]';
% b_RLS = theta_RLS(end/2+1 : end)';
% G_RLS = tf(b_RLS, a_RLS, sampling_time_02)
%  % Simulate the output of the system using the input data 
%  % and the transfer function G
% time = (0 : length(input_02) - 1)' * sampling_time_02;
% simulated_output = lsim(G_RLS, input_02, time);
% 
% 
% figure
% plot(time, simulated_output)
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to Training signal, performed in the operating point, Ts = 0.4')
% legend('RLS model response', 'Plant response')
% hold on;
% plot(time_training,output_02)
% legend('RLS model response', 'Plant response')
% 
% 
% 
% validation_Model_output_RLS_02 = lsim(G_RLS, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
% figure
% plot(validation_Model_time, validation_Model_output_RLS_02)
% % Ploting output of simulation
% hold on;
% % plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to validation signal, performed in the operating point, Ts = 0.4 s')
% legend('RLS model response', 'Plant response')
% 
% 
% % Calculating mean of abs error
% d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_RLS_02));
% fprintf("RLS: MAE = %.3f; Ts = 0.4 s \n",d)
% 
% %% STOHASTIC APROXIMATION
% theta_SA = stochastic_approximation(input_02, output_02, order, delay);
% a_SA = [1; theta_SA(1 : order)]';
% b_SA = theta_SA(order + 1 : 2*order)';
% G_SA = tf(b_SA, a_SA, sampling_time_02)
% 
%  % Simulate the output of the system using the input data 
% simulated_output_SA = lsim(G_SA, input_02, time);
% 
% figure
% plot(time, simulated_output_SA)
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to Training signal, performed in the operating point, Ts = 0.4')
% legend('SA model response', 'Plant response')
% hold on;
% plot(time_training,output_02)
% legend('SA model response', 'Plant response')
% 
% validation_Model_output_SA_02 = lsim(G_SA, validation_input(signal_duration1/sampling_time_02:end) - amplitude_min, validation_time(signal_duration1/sampling_time_02:end)) + output_02_oo;
% figure
% plot(validation_Model_time, validation_Model_output_SA_02)
% % Ploting output of simulation
% hold on;
% % plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end) - simulated_validation_output_oo)
% plot(validation_Model_time,simulated_validation_output(signal_duration1/sampling_time_02:end))
% grid on;
% xlabel('t / s')
% ylabel('U / V')
% title('Response of the plant to validation signal, performed in the operating point, Ts = 0.4')
% legend('SA model response', 'Plant response')
% 
% 
% % Calculating sum of error
% d = mean(abs(simulated_validation_output(signal_duration1/sampling_time_02:end) - validation_Model_output_SA_02));
% fprintf("SA: MAE = %.3f; Ts = 0.4 s \n",d)
% 
% 
% 
% %% ETFE
% 
% 
% 
% % Number of samples
% N_02 = length(output_02);
% 
% % Frequency resolution
% F_02 = 1 / (N_02 * sampling_time_02 );
% 
% % Maximum frequency (1/2 of sampling frequency)
% f_max_02 = 1 / (2 * sampling_time_02);
% 
% % Calculte the FT of input and output signals
% input_02_fft = fft(input_02)* sampling_time_02;
% output_02_fft = fft(output_02') * sampling_time_02;
% 
% % Calculate the ETFE
% G_bode_02 = output_02_fft ./ input_02_fft;
% G_bode_02_half = G_bode_02(1 : N_02/2);
% 
% 
% % Frequencies used for Bode plot
% frequencies_02 = 0 : F_02 : f_max_02 - F_02;
% 
% % Get frequency response of the model
% [amplitude_02, phase_02] = bode(G_02, 2 * pi * frequencies_02);
% amplitude_02 = amplitude_02(:);
% phase_02 = phase_02(:);
% 
% 
% 
% G_bode_02_half_phase = angle(G_bode_02_half);
% % G_bode_02_half_phase_unwrapped = unwrap(G_bode_02_half_phase);
% G_bode_02_half_phase_unwrapped = custom_unwrap(G_bode_02_half_phase,7,pi);
% 
% 
% %% SPECTRAL ANALYSIS
% 
% 
% phi_uu = my_cross_correlation(input_02', input_02');
% phi_uy = my_cross_correlation(input_02', output_02');
% T_m = length(phi_uu)/8;
% w = custom_hann(length(phi_uu), 2*T_m);
% 
% Phi_uu = fft(phi_uu .* w);
% Phi_uy = fft(phi_uy .* w);
% 
% 
% G_spectral_analysis = Phi_uy'./Phi_uu';
% 
% 
% 
% % Number of samples
% N = length(G_spectral_analysis);
% 
% % Frequency resolution
% F = 1 / (N * sampling_time_02 );
% 
% % Maximum frequency (1/2 of sampling frequency)
% f_max = 1 / (2 * sampling_time_02);
% 
% % Frequencies used for Bode plot
% frequencies = 0 : F : f_max - F;
% 
% 
% % Bode plot
% figure
% subplot(2,1,1)
% set(gca, 'XScale', 'log')
% hold on
% grid on
% SA_plot = semilogx(frequencies, 20*log10(abs(G_spectral_analysis(1:end/2))));
% LS_plot = semilogx(frequencies_02 , 20 * log10(amplitude_02),'g' );
% ETFE = semilogx(frequencies_02 , 20 * log10(abs(G_bode_02_half)),'r:');
% for i = 1 : length(harmonic_frequencies)
%     scatter_handles(i) = scatter(harmonic_frequencies(i), 20 * log10(frequency_response(i, 1)), 'filled', 'k');
% end
% ylabel('|G(j\omega)| / dB')
% legend([SA_plot,LS_plot,ETFE, scatter_handles(1)], 'Spectral analysis', 'LS','ETFE', 'Orthogonal Correlation');
% title('Amplitude response (Ts = 0.4)')
% 
% G_spectral_analysis_phase = angle(G_spectral_analysis);
% G_spectral_analysis_phase_unwraped = custom_unwrap(G_spectral_analysis_phase,3,pi);
% %G_spectral_analysis_phase_unwraped(2) = G_spectral_analysis_phase_unwraped(2) + 2*pi;
% % G_bode_02_half_phase_unwrapped(3:7) = G_bode_02_half_phase_unwrapped(3:7) +pi
% % G_bode_02_half_phase_unwrapped(7) = G_bode_02_half_phase_unwrapped(7) + 2*pi
% 
% subplot(2,1,2)
% set(gca, 'XScale', 'log')
% hold on
% grid on
% semilogx(frequencies , rad2deg(G_spectral_analysis_phase_unwraped(1:end/2)));
% hold on;
% semilogx(frequencies_02 , (phase_02),'g')
% hold on;
% % G_bode_02_half_phase_unwrapped(13) = G_bode_02_half_phase_unwrapped(13) - pi;
% semilogx(frequencies_02 , rad2deg(G_bode_02_half_phase_unwrapped),'r:')
% hold on
% for i = 1 : length(harmonic_frequencies)
%     scatter(harmonic_frequencies(i), frequency_response(i, 2), 'filled', 'k')
% end
% xlabel('f / Hz')
% ylabel('\theta(G(j\omega)) / deg')
% title('Phase response')

