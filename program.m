%Caner Budak - 280070846121 - el22606

%% for question 1.1, i) 
%built in fucntion fir1


% Moving average filter with N = 2
N1 = 2;
b1 = fir1(N1, 1/N1, 'low');
a1 = 1;
b = fir1(N, Wn, 'type');
% Moving average filter with N = 4
N2 = 4;
b2 = fir1(N2, 1/N2, 'low');
a2 = 1;

% Moving average filter with N = 10
N3 = 10;
b3 = fir1(N3, 1/N3, 'low');
a3 = 1;

% Generate a test signal
x = randn(1, 100);

% Filter the signal using the moving average filter with N = 2
y1 = filter(b1, a1, x);

% Filter the signal using the moving average filter with N = 4
y2 = filter(b2, a2, x);

% Filter the signal using the moving average filter with N = 10
y3 = filter(b3, a3, x);

%% for question 1.1 ii) 
% Moving average filter with N = 2
N1 = 2;
b1 = fir1(N1, 1/N1, 'low');
a1 = 1;

% Moving average filter with N = 4
N2 = 4;
b2 = fir1(N2, 1/N2, 'low');
a2 = 1;

% Moving average filter with N = 10
N3 = 10;
b3 = fir1(N3, 1/N3, 'low');
a3 = 1;

% Plot the amplitude and phase response of the filters
[h1, w1] = freqz(b1, a1);
[h2, w2] = freqz(b2, a2);
[h3, w3] = freqz(b3, a3);

figure;
subplot(2,1,1);
plot(w1/pi, abs(h1), 'b', w2/pi, abs(h2), 'g', w3/pi, abs(h3), 'r');
legend('N = 2', 'N = 4', 'N = 10');
title('Magnitude Response');
xlabel('Normalized Frequency (x pi radians/sample)');
ylabel('Magnitude');

subplot(2,1,2);
plot(w1/pi, angle(h1), 'b', w2/pi, angle(h2), 'g', w3/pi, angle(h3), 'r');
legend('N = 2', 'N = 4', 'N = 10');
title('Phase Response');
xlabel('Normalized Frequency (x pi radians/sample)');
ylabel('Phase (radians)');
%% for question 1.1 iii) 
% Moving average filter with N = 2
N1 = 2;
b1 = fir1(N1, 1/N1, 'low');
a1 = 1;

% Moving average filter with N = 4
N2 = 4;
b2 = fir1(N2, 1/N2, 'low');
a2 = 1;

% Moving average filter with N = 10
N3 = 10;
b3 = fir1(N3, 1/N3, 'low');
a3 = 1;

% Convert the filters to zero-pole representation
[z1, p1, k1] = tf2zp(b1, a1);
[z2, p2, k2] = tf2zp(b2, a2);
[z3, p3, k3] = tf2zp(b3, a3);

% Plot the pole-zero diagram of the filters
figure;
subplot(3,1,1);
zplane(z1, p1);
title('Pole-Zero Diagram (N = 2)');

subplot(3,1,2);
zplane(z2, p2);
title('Pole-Zero Diagram (N = 4)');

subplot(3,1,3);
zplane(z3, p3);
title('Pole-Zero Diagram (N = 10)');

%% for question 1.2 i) 
%also I would like to mention in here, I couldn't find "zpd.m" function in
%supplementary materials 

% Step 1: Define the poles and zeros
p = [0.68+0.51i 0.68-0.51i];
z = [1.2 ; -0.6];

% Step 2: Calculate the transfer function coefficients
[b,a] = zp2tf(z,p,1);

% Step 3: Plot the system's pole-zero diagram
zplane(b,a)

%% for question 1.2 ii) 
% Compute the frequency response
[h,w] = freqz(b,a);

% Plot the amplitude response
figure;
plot(w, abs(h));
title('Amplitude response');
xlabel('Normalized frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Plot the phase response
figure;
plot(w, angle(h));
title('Phase response');
xlabel('Normalized frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

%% for question 1.2 iii) 

% Compute the impulse response
imp = impz(b,a);

% Compute the step response
step = stepz(b,a);

% Plot the impulse response
figure;
stem(imp);
title('Impulse response');
xlabel('Sample index');
ylabel('Amplitude');

% Plot the step response
figure;
stem(step);
title('Step response');
xlabel('Sample index');
ylabel('Amplitude');

%% for question 1.2 iv)

% Original filter
z = [1.2, -0.6];
p = [0.68+0.51i, 0.68-0.51i];
k = 1;
sys = zpk(z, p, k);
[b, a] = zp2tf(z, p, k);

% Plot step response
step = stepz(b, a);
figure;
stem(step);
title('Step response (original filter)');
xlabel('Sample index');
ylabel('Amplitude');

% Plot amplitude response
[h, w] = freqz(b, a);
figure;
plot(w, abs(h));
title('Amplitude response (original filter)');
xlabel('Normalized frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% New pole locations
p1 = [0.76+0.57i, 0.76-0.57i];
p2 = [0.8+0.6i, 0.8-0.6i];
p3 = [0.84+0.63i, 0.84-0.63i];

% Move poles to new locations and plot step response
sys1 = zpk(z, p1, k);
[b1, a1] = zp2tf(z, p1, k);
step1 = stepz(b1, a1);
figure;
stem(step1);
title('Step response (pole locations: 0.76±0.57i)');
xlabel('Sample index');
ylabel('Amplitude');

% Move poles to new locations and plot amplitude response
sys1 = zpk(z, p1, k);
[b1, a1] = zp2tf(z, p1, k);
[h1, w1] = freqz(b1, a1);
figure;
plot(w1, abs(h1));
title('Amplitude response (pole locations: 0.76±0.57i)');
xlabel('Normalized frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Move poles to new locations and plot step response
sys2 = zpk(z, p2, k);
[b2, a2] = zp2tf(z, p2, k);
step2 = stepz(b2, a2);
figure;
stem(step2);
title('Step response (pole locations: 0.8±0.6i)');
xlabel('Sample index');
ylabel('Amplitude');

% Move poles to new locations and plot step response
sys3 = zpk(z, p3, k);
[b3, a3] = zp2tf(z, p3, k);
step3 = stepz(b3, a3);
figure;
stem(step3);
title('Step response (pole locations: 0.84±0.63i)');
xlabel('Sample index');
ylabel('Amplitude');

%% for question 1.2 v)
% Original filter
z = [1.2, -0.6];
p = [0.68+0.51i, 0.68-0.51i];
k = 1;
sys = zpk(z, p, k);

% Generate input pulse with T=50s and simulate system response
T1 = 50;
t = 0:1:100;
pulse1 = gensig('pulse', T1, 100, 1);
response1 = lsim(sys, pulse1, t);

% Plot input pulse and system response
figure;
subplot(2,1,1);
plot(t, pulse1);
title('Input pulse (T=50s)');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, response1);
title('System response (T=50s)');
xlabel('Time (s)');
ylabel('Amplitude');

% Generate input pulse with T=5s and simulate system response
T2 = 5;
pulse2 = gensig('pulse', T2, 100, 1);
response2 = lsim(sys, pulse2, t);

% Plot input pulse and system response
figure;
subplot(2,1,1);
plot(t, pulse2);
title('Input pulse (T=5s)');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, response2);
title('System response (T=5s)');
xlabel('Time (s)');
ylabel('Amplitude');

%% for question 1.2 vi) 
z = [1.2, -0.6];
p = [0.8i, -0.8i];
k = 1;
sys = zpk(z, p, k);
[b, a] = zp2tf(z, p, k);
freqz(b, a);

% Compute the frequency response
[h,w] = freqz(b,a);

% Plot the amplitude response
figure;
plot(w, abs(h));
title('Amplitude response');
xlabel('Normalized frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Plot the phase response
figure;
plot(w, angle(h));
title('Phase response');
xlabel('Normalized frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

%% for question 2.1 i) 

[flute_note, fs] = audioread('flute_note.wav');
sound(flute_note, fs);

[clarinet_note, fs] = audioread('clarinet_note.wav');
sound(clarinet_note, fs);

[cello_note, fs] = audioread('cello_note.wav');
sound(cello_note, fs);

%% for question 2.1 ii) 
% Flute note
t_flute = (0:length(flute_note)-1)/fs;
hold on 
plot(t_flute(1:22050), flute_note(1:22050));
xlabel('Time (s)');
ylabel('Amplitude');
title('Flute Note');

% Clarinet note
t_clarinet = (0:length(clarinet_note)-1)/fs;
plot(t_clarinet(1:22050), clarinet_note(1:22050));
xlabel('Time (s)');
ylabel('Amplitude');
title('Clarinet Note');

% Cello note
t_cello = (0:length(cello_note)-1)/fs;
plot(t_cello(1:22050), cello_note(1:22050));
xlabel('Time (s)');
ylabel('Amplitude');
title('Cello Note');
hold off

%% for question 2.1 iii) 
% Flute note
N = length(flute_note);
f = (0:N-1)*(fs/N);
F = fft(flute_note, N);
figure;
hold on
plot(f, abs(F));
xlim([0 fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Flute Note Spectrum');

% Clarinet note
N = length(clarinet_note);
f = (0:N-1)*(fs/N);
F = fft(clarinet_note, N);
figure;
plot(f, abs(F));
xlim([0 fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Clarinet Note Spectrum');

% Cello note
N = length(cello_note);
f = (0:N-1)*(fs/N);
F = fft(cello_note, N);
figure;
plot(f, abs(F));
xlim([0 fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Cello Note Spectrum');
hold off

%% for question 2.1 iv) 
% Normalize the signals
flute_note = flute_note / max(abs(flute_note));
clarinet_note = clarinet_note / max(abs(clarinet_note));
cello_note = cello_note / max(abs(cello_note));

% Define the window size
Nwin = 1000;
w = ones(Nwin, 1);

% Calculate the energy using a sliding window with conv() command
E_flute = conv(abs(flute_note).^2, w, 'same');
E_clarinet = conv(abs(clarinet_note).^2, w, 'same');
E_cello = conv(abs(cello_note).^2, w, 'same');

% Plot the results
t = (0:length(flute_note)-1)/fs;
figure;
plot(t, flute_note, t, E_flute);
xlabel('Time (s)');
ylabel('Amplitude/Energy');
title('Flute Note and Energy');
legend('Flute Note', 'Energy');

t = (0:length(clarinet_note)-1)/fs;
figure;
plot(t, clarinet_note, t, E_clarinet);
xlabel('Time (s)');
ylabel('Amplitude/Energy');
title('Clarinet Note and Energy');
legend('Clarinet Note', 'Energy');

t = (0:length(cello_note)-1)/fs;
figure;
plot(t, cello_note, t, E_cello);
xlabel('Time (s)');
ylabel('Amplitude/Energy');
title('Cello Note and Energy');
legend('Cello Note', 'Energy');

%% for question 2.1 v) 

% Load and listen the signal
[cello_note_noisy, fs] = audioread('cello_note_noisy.wav');
sound(cello_note_noisy, fs);

% Plot the spectrum
N = length(cello_note_noisy);
X = fft(cello_note_noisy);
freq = (0:N-1)/N*fs;
figure;
plot(freq, abs(X));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Noisy Cello Note');

%% for question 2.1 vi) 

% Define filter length
filter_length = 50;

% Create filter coefficients for a moving average filter
b = ones(1, filter_length)/filter_length;

% Apply the filter to the noisy signal
cello_note_filtered = filter(b, 1, cello_note_noisy);

% Plot the spectrum of the filtered signal
N = length(cello_note_filtered);
X = fft(cello_note_filtered);
freq = (0:N-1)/N*fs;
figure;
plot(freq, abs(X));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Filtered Cello Note');

% Listen to the filtered signal
sound(cello_note_filtered, fs);

%% for question 2.1 vii) 
% Define filter parameters
fc = 9*f0;   % center frequency of the passband
bw = 50;     % bandwidth of the passband

% Create filter coefficients for a band-pass filter
[b,a] = butter(4, [2*fc/fs-bw/2, 2*fc/fs+bw/2], 'bandpass');

% Apply the filter to the cello note
cello_note_filtered_9 = filter(b, a, cello_note);

% Plot the spectrum of the filtered signal
N = length(cello_note_filtered_9);
X = fft(cello_note_filtered_9);
freq = (0:N-1)/N*fs;
figure;
plot(freq, abs(X));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Cello Note, 9th Harmonic Filtered');

% Plot a segment of the filtered signal
t = (0:N-1)/fs;
figure;
plot(t, cello_note_filtered_9);
xlabel('Time (s)');
ylabel('Amplitude');
title('Cello Note, 9th Harmonic Filtered');

%% for question 2.2 i) 
% I prefere to use cello_note signal for segmentation
% Load the cello note signal
[x, fs] = audioread('cello_note.wav');

% Calculate the spectrum of the signal
N = length(x);
X = fft(x);
f = (0:N-1)*(fs/N);

% Find the peak frequency
[~, idx] = max(abs(X));
fundamental_freq = f(idx);
period = 1/fundamental_freq;

% Select a segment of length 15 times the period
segment_length = round(15*period);
start_index = round(N/2 - segment_length/2);
end_index = start_index + segment_length - 1;
segment = x(start_index:end_index);

% Plot the segment and the original signal
t = (0:N-1)/fs;
plot(t, x, t(start_index:end_index), segment);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Selected Segment');

%% for question 2.2 ii) 
% Segment can be input variable from previous question
% Calculate the Fourier transform of the selected segment
X = fft(segment);
N = length(segment);
f = (0:N-1)*(fs/N);

% Plot the amplitude of the Fourier transform
figure;
plot(f, abs(X)/N);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude of Fourier Transform');

% Calculate the amplitude coefficients for each harmonic
threshold = 0.1; % set an amplitude threshold for significant harmonics
a1 = max(abs(X));
c = abs(X)/a1;
c(c < threshold) = 0; % zero out coefficients below threshold

%% for question 2.2 iii) 
% We can use the phase vector c from previous question 
% Find the indices of significant harmonics
harmonics = find(c);

% Calculate the phase angles for each harmonic
phi = angle(X(harmonics));


%% for question 2.2 iv) 

% Create a time vector with fs samples per second
t = linspace(0, N/fs, N);

% Initialize the reconstructed signal
x_recon = zeros(size(t));

% Add up the individual sinusoids
for i = 1:length(harmonics)
    x_recon = x_recon + c(harmonics(i)) * cos(2*pi*harmonics(i)*t + phi(i));
end

%% for question 2.2 v) 

% Select a segment of the signal with M samples, where M is a multiple of the fundamental period
M = round(T / dt) * 15;
segment = x(1:M);

% Compute the DFT of the segment
X = fft(segment);

% Find the significant harmonics
threshold = 0.2;  % adjust as necessary
harmonics = find(abs(X) >= threshold * max(abs(X)));
c = abs(X(harmonics)) / abs(X(1));
phi = angle(X(harmonics));

% Reconstruct the signal using the significant harmonics
t = linspace(0, M*dt, M);
x_recon = zeros(size(t));
for i = 1:length(harmonics)
    x_recon = x_recon + c(i) * cos(2*pi*harmonics(i)*t + phi(i));
end

% Normalize the signals to [-1, 1]
segment_norm = segment / max(abs(segment));
x_recon_norm = x_recon / max(abs(x_recon));

% Compare the original and reconstructed signals
figure;
plot(t, segment_norm, 'b-', t, x_recon_norm, 'r--');
legend('Original', 'Reconstructed');
xlabel('Time (s)');
ylabel('Amplitude');

% Listen to the reconstructed signal
x_recon_norm = x_recon_norm / max(abs(x_recon_norm));  % normalize to [-1, 1]
sound(x_recon_norm, fs);

% Compare the results using different numbers of harmonics
figure;
subplot(3, 1, 1);
plot(t, segment_norm, 'b-', t, x_recon_norm, 'r--');
legend('Original', 'Reconstructed');
xlabel('Time (s)');
ylabel('Amplitude');
title('All harmonics');

harmonics = 1;
c = abs(X(harmonics)) / abs(X(1));
phi = angle(X(harmonics));
x_recon = zeros(size(t));
for i = 1:length(harmonics)
    x_recon = x_recon + c(i) * cos(2*pi*harmonics(i)*t + phi(i));
end
x_recon_norm = x_recon / max(abs(x_recon));
subplot(3, 1, 2);
plot(t, segment_norm, 'b-', t, x_recon_norm, 'r--');
legend('Original', 'Reconstructed');
xlabel('Time (s)');
ylabel('Amplitude');
title('Fundamental harmonic only');

harmonics = 1:5;
c = abs(X(harmonics)) / abs(X(1));
phi = angle(X(harmonics));
x_recon = zeros(size(t));
for i = 1: length(harmonics)
    x_recon = x_recon + c(i) * cos(2*pi*harmonics(i)*t + phi(i));
end

%% for question 2.2 vi) 

fs = 44100; % Sampling frequency to be like original one 
filename = 'reconstructed_cello_note.wav';
audiowrite(filename, reconstructed_signal, fs);


