%% Non-GUI ALEX data processing 
% Got tired of using smViewer for data processing so wrote a function to do
% it automatically here. A couple assumptions
% 
% 1. Works for data that is collected under ALEX 
% 2. Assumes channel 2 is RR and channel 1 is GG (default now in
% MicroManager). 
% 3. Maps files directly and does not use an alignment file. 

% David S. White 
% 2023-05-04
% -------------------------------------------------------------------------
% load data