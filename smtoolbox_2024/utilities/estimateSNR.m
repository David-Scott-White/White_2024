function [snr, snb, sigma] = estimateSNR(timeSeries)
% Estimate signal to noise ratio (SNR) of a trace
% sigma is estimated via estimateNoise.m 
% Signal is detected via kmeans clustering of two states

% SNR = (signal intensity state 2)/sigma
% SNB = (signal intensity state 2) - (signal intensity state 1)/sigma
% sigma = value from estimateNoise 

% input: timeSeries

% David S. White 
% 2022-02-10
% MIT 
% -------------------------------------------------------------------------


states = kmeansElkan(timeSeries,2); 
sigma = estimateNoise(timeSeries); 
snr = max(states)/sigma; 
snb = abs(diff(states))/sigma; 