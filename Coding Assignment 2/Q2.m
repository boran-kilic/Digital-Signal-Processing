clear; close all; clc;
%% Downsampling

[audio, fs_orig] = audioread('original_audio.wav');

fs_target = 20000;
[p, q] = rat(fs_target/fs_orig);  % 20/48 = 5/12


% Up-sample
audio_up = upsample(audio, p); 

% Anti-aliasing filter before downsampling
filt = fir1(100, 1/q); 
audio_filt = conv(audio_up, filt, 'same');

% Down-sample
audio_20k = downsample(audio_filt, q); 
audiowrite('audio_20kHz.wav', audio_20k, fs_target); 

%% FIR Bandpass Filter Windowing
        
L = 10;           
n = -(L-1)/2 : (L-1)/2;   
n1 = 0:L-1;

f_low = 275;         
f_high = 3250;       

w_low = 2*pi*f_low/fs_target;
w_high = 2*pi*f_high/fs_target;

h_ideal = (sin(w_high*n) - sin(w_low*n)) ./ (pi*n);
h = h_ideal;

%% Filtering
audio_filtered = conv(audio_20k, h, 'same');
audiowrite('filtered_audio_20kHz.wav', audio_filtered, fs_target); 

%% Number of Samples
fprintf("Number of Samples\n");
fprintf('Original recording: %d\n', length(audio));           % 480,000
fprintf('Resampled recording: %d\n', length(audio_20k)); % 200,000
fprintf('Recording after filtering: %d\n', length(audio_filtered)); % 200,000
%% Frequency Domain Analysis
Fs1 = 48000;
n1 = length(audio);
f1 = linspace(-Fs1/2, Fs1/2, n1);
Y1 = fftshift(abs(fft(audio)));

figure;
plot(f1, Y1);
title('FT of Original Audio (48 kHz)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;

Fs2 = 20000;
n2 = length(audio_20k);
f2 = linspace(-Fs2/2, Fs2/2, n2);
Y2 = fftshift(abs(fft(audio_20k)));

figure;
plot(f2, Y2);
title('FT of Resampled Audio (20 kHz)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;

Fs3 = 20000;  
n3 = length(audio_filtered);
f3 = linspace(-Fs3/2, Fs3/2, n3);
Y3 = fftshift(abs(fft(audio_filtered)));

figure;
plot(f3, Y3);
title('FT of Filtered Audio');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;
