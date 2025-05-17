clear; close all; clc;
%% Audio Recording
fs_orig = 48000;
recObj = audiorecorder(fs_orig, 16, 1);
duration = 11;

disp('Start speaking...');
recordblocking(recObj, duration); 
disp('End of Recording.');

audio = getaudiodata(recObj);
audio = audio(fs_orig-1:fs_orig*11); 
audiowrite('original_audio.wav', audio, fs_orig);