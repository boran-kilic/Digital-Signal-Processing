clear; close all; clc;
%%
% Parameters
fs = 20000;              
L = 10;                  
f_low =300;         
f_high = 3000;  
f = [0 0 f_low-299 f_high f_high+100 fs/2]/ (fs/2);     
a = [0 0 1 1 0 0];                  

h = firpm(L-1, f, a);
[H_temp, ~] = freqz(h, 1, 1024);
h = h / max(abs(H_temp));  
[H, freq] = freqz(h, 1, 1024, fs);

% Impulse Response
figure;
stem(0:L-1, h, 'filled');
xlabel('n'); ylabel('h[n]');
title('Impulse Response');
grid on;

% Magnitude Response
figure;
plot(freq, abs(H));
xlabel('Frequency (Hz)');
ylabel('|H(f)|');
title('Magnitude Response (Minimax)');
grid on; xlim([0 fs/2]);
hold on;
yline(0.707, 'r--', '-3 dB');
hold off 

% Phase Response
figure;
plot(freq, unwrap(angle(H)));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Response');
grid on;

% Pole-Zero Plot
figure;
zplane(roots(h),[]);
title('Pole-Zero Plot');

% Z-Transform
syms z
H_z = 0;
for k = 1:length(h)
    H_z = H_z + h(k) * z^(-(k-1));
end

