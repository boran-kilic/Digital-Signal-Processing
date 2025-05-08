close all
clear 
clc

rng(2,"twister")
N = 9;
real_part = randn(1,N);
imag_part = randn(1,N);
x = real_part + 1i*imag_part;

fprintf('N = %d \n', N);
n = 0:N-1; 

figure;
subplot(2,1,1);
stem(n, real(x), 'filled');
title('Real Part of x[n]');
xlabel('n');
ylabel('Re\{x[n]\}');
xlim([0 N-1]);
grid on;

subplot(2,1,2);
stem(n, imag(x), 'filled');
title('Imaginary Part of x[n]');
xlabel('n');
ylabel('Im\{x[n]\}');
xlim([0 N-1]);
grid on;

t1 = tic;
Xsum = DFT_summation(x);
elapsed_time = toc(t1);
fprintf('Elapsed time for DFT_summation is %.6f seconds.\n', elapsed_time);
plotDFT(Xsum, 'Summation Formula')

t2 = tic;
X_9pt_fft = nineptFFT(x);
elapsed_time = toc(t2);
fprintf('Elapsed time for 9-pt FFT is %.6f seconds.\n', elapsed_time);
plotDFT(X_9pt_fft, '9-pt FFT')

t3 = tic;
X_fft = fft(x);
elapsed_time = toc(t3);
fprintf('Elapsed time for fft is %.6f seconds.\n', elapsed_time);
plotDFT(X_fft, 'FFT')

fprintf('\n');
disp(['Summation Formula:   ', num2str(norm(X_fft - Xsum))]);
disp(['DFT matrix:          ', num2str(norm(X_fft - X_9pt_fft))]);
disp(['fft:                 ', num2str(norm(X_fft - X_fft))]);
fprintf('\n');

function X = nineptFFT(x)
N = length(x);
X = zeros(1, N);
for r = 0:(N/3 - 1)
    for p = 0:2
        temp = 0;
        for n = 0:(N/3 - 1)
            inner_sum = 0;
            for l = 0:2
                inner_sum = inner_sum + x(n + l*N/3 + 1) * exp(-1j*2*pi*p*l/3);
            end
            temp = temp + inner_sum * exp(-1j*2*pi*n*r/(N/3)) * exp(-1j*2*pi*p*n/N);
        end
        X(3*r + p + 1) = temp;
    end
end
end


function X = DFT_summation(x)
N = length(x);
X = zeros(1, N); 
for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1) * exp(-1i * 2 * pi * k * n / N);
    end
end
end


function plotDFT(X, str)
N = length(X);
n = 0:N-1;

figure;
subplot(2,1,1);
stem(n, abs(X), 'filled');
title(['Magnitude of X[k] using ', str]);
xlabel('k');
ylabel('|X[k]|');
xlim([0 N-1]);
grid on;

subplot(2,1,2);
stem(n, angle(X), 'filled');
title(['Phase of X[k] using ', str]);
xlabel('k');
ylabel('∠X[k] (radians)');
xlim([0 N-1]);
grid on;
end


