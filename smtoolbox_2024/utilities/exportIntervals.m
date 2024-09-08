function export = exportIntervals(Intervals, targetSoftware)
% -------------------------------------------------------------------------
% Convert Intervals.mat from imscroll to data structures for other programs
% -------------------------------------------------------------------------
%
% Input:
%   intervals = "Intervals.mat" file from imscroll.
%       > if empty, uigetfile
%   targetSoftware  = intended software to read intervals. Currently supported:
%       > None = returned as simple data structure.
%       > DISCO (White et al. 2020)
%       > ebFRET (van de Meent et al., 2013) (also vbFRET)
%           *** Saves as stacked donr-acceptor traces dat format ***
%
% Output:
%   export. Variable output depending on the targetSoftware selected.
%       "None" = struct;
%       "DISCO" returns export.data struct where data = DISCO struct;
%       "ebFRET" returns export an array [N, 4].
%
% Dependencies: 
%   computeCenters.m (DISCO only, in DISCO package)
%
% Author: David S. White
% Updated: 2021-10-06
% License: MIT
% -------------------------------------------------------------------------
export = [];
if ~exist('Intervals', 'var') || isempty(Intervals)
    [loadFile, loadPath] = uigetfile('*.mat', 'Load and Interval file','Multiselect', 'off');
    temp = load([loadPath, loadFile]);
    if isfield(temp, 'Intervals')
        Intervals = temp.Intervals;
    else
        disp('Error. imscroll Intervals.mat file not loaded.');
        return
    end
end
if ~exist('targetSoftware', 'var') || isempty(targetSoftware)
    targetSoftware = 'ebFRET';
end

traces = Intervals.AllTracesCellArray;
nTraces = size(traces,1);
export = struct;
for i = 1:nTraces
    trace1 = traces(i,12);
    trace2 = traces(i,9);
    time_s = trace2{1};
    
    frames = trace1{1}(:,1);
    intensity = trace1{1}(:,2);
    
    % Correct for values below zero...
    minInt = min(intensity);
    if min(minInt) < 0
        intensity = intensity + abs(minInt);
    end
    
    ideal = Intervals.AllTracesCellArray(i, 13);
    ideal = ideal{1}(:,1);
    
    % For some reason the entire last event if dropped....
    while length(ideal) < length(intensity)
        ideal = [ideal; ideal(end)];
    end
    
    while length(time_s) > length(intensity)
        time_s(end) = []; 
    end
    
    switch targetSoftware
        case 'None'
            export(i,1).frames = frames;
            export(i,1).time_s = time_s;
            export(i,1).intensity = intensity;
            export(i,1).ideal = ideal;
            
        case 'DISCO'
            ideal(ideal==-3) = 1; % first event bound
            ideal(ideal== 3) = 1; % last event bound
            ideal(ideal==-2) = 0; % first event unbound
            ideal(ideal== 2) = 0; % last event unbound
            
            [components,ideal, class] = computeCenters(intensity, ideal+1);
            
            export.data.rois(i,1).time_series = intensity;
            export.data.rois(i,1).frames = frames;
            export.data.rois(i,1).time_s = time_s;
            export.data.rois(i,1).disc_fit.components = components;
            export.data.rois(i,1).disc_fit.ideal = ideal;
            export.data.rois(i,1).disc_fit.class = class;
            
        case 'ebFRET' % assume "intensity" is the FRET trace and make two fake traces
            export(i,1).frames = frames;
            export(i,1).time_s = time_s;
            export(i,1).intensity = intensity;
            
            % E = A/(A+D); Assume intensity = Acceptor (A), therefore
            % create a scaled Donor (D) by D = max(A) - A.
            % See Methods of (Goldschen-Ohm eLife 2016)
            
            % min-max normalize A to scale intensities for FRET
            A = intensity;
            A = (A-min(A))./(max(A)-min(A));
            
            A = A + 1e-3; % non-zero values needed
            A = A * 100;  % add some noise back to it. 
            D = max(A)-A;
            E = A./(A+D);
            export(i,1).acceptor = A;
            export(i,1).donor = D;
            export(i,1).fret = E;
            
    end
end

if exist('loadPath', 'var')
    tempName = [loadPath, loadFile(1:end-4), '_', targetSoftware, '.mat'];
    [saveFile, savePath] = uiputfile(tempName);
else
    [saveFile, savePath] = uiputfile('*.mat');
end

if saveFile
    switch targetSoftware
        case 'ebFRET'
            m = length(export);
            n = length(export(1).donor);
            traces = zeros(n+1, m*2);
            idx = 1;
            for i = 1:2:2*m
                traces(:, i) = [i;export(idx).donor];
                traces(:, i+1) = [i;export(idx).acceptor];
                idx = idx+1;
            end
            writematrix(traces, [savePath, saveFile(1:end-4), '.dat']);
            disp([saveFile(1:end-4), '.dat saved.'])
            
        case 'DISCO'
            data = export.data;
            save([savePath, saveFile], 'data');
            disp([saveFile, ' saved.'])
            
        case 'None'
            save([savePath, saveFile], 'export');
            disp([saveFile, ' saved.'])
            
    end
end

