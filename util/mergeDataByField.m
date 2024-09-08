function analysis_out  = mergeDataByField(analysis_in, field_name, varargin)
% David S. White 
% 2023-06-14

% works on .mat structures generated from createDataset
bootstrap = 0; 
useFitLimits = []; 
dwellLimits = [];
dropLastUnbound = 0; 

switch field_name
    case{'rna_nM'}
        field_values = vertcat(analysis_in.rna_nM);
    case{'U1C_nM', 'U1C'}
        field_values = vertcat(analysis_in.U1C_nM);
    case {'cpd_uM'}
        field_values = vertcat(analysis_in.cpd_uM);
end

for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'bootstrap'
            bootstrap = varargin{i+1};
        case 'fitLimits'
            useFitLimits = varargin{i+1};
        case 'tauGuess'
            tauGuess = varargin{i+1};
        case 'dropLastUnbound'
            dropLastUnbound = varargin{i+1};
    end
end

% just merge dwells, fraction bound, and refit? 

[unique_field_values] = unique(field_values);
for i = 1:length(unique_field_values)
    idx = find(field_values == unique_field_values(i));

    a = struct; 
    a.field_value = []; 
    a.field_name = [];
    a.nfiles = 0; 
    a.nrois = 0;
    a.frameRate_s = {}; 
    a.nFrames = {}; 
    a.time_s = {}; 
    a.events = {}; 
    a.nevents = 0; 
    a.timeseries = {}; 
    a.stateseq = {}; 
    a.dwells = cell(2,1); 
    a.dwellsCDF = cell(2,1); 
    a.dwellsmle = cell(2,1); 
    a.timeToFirst = []; 
    a.timeToFirstMLE = [];
    a.timetofirstmle = [];
    a.timetofirstmle_se = [];
    a.unboundmle = []; 
    a.unboundmle_se = []; 
    a.unbound_k = 0; 
    a.boundmle = [];
    a.boundmle_se = []; 
    a.bound_k = 0; 
    a.fractionbound = []; 
    a.fractionbound_fit = []; 
    
    for j = 1:numel(idx)
        a.nfiles = a.nfiles + 1;
        a.nrois = a.nrois + analysis_in(idx(j)).nrois;
        
        a.frameRate_s{j} = analysis_in(idx(j)).frameRate_s;
        a.nFrames{j} = analysis_in(idx(j)).nframes;
        a.time_s{j} = analysis_in(idx(j)).time_s;
        a.timeseries{j} = analysis_in(idx(j)).timeseries;
        a.stateseq{j} = analysis_in(idx(j)).stateseq;
        if dropLastUnbound
            for m = 1:length(analysis_in(idx(j)).events)
                if size(analysis_in(idx(j)).events{m},1) > 1
                    if analysis_in(idx(j)).events{m}(end,4) == 0
                        dt =  analysis_in(idx(j)).events{m}(end,3)*analysis_in(idx(j)).frameRate_s;
                        analysis_in(idx(j)).events{m}(end,:) = [];
                        analysis_in(idx(j)).nevents = analysis_in(idx(j)).nevents - 1;
                        bidx = analysis_in(idx(j)).events{m}(:,4) > 0;
                        bf = sum(analysis_in(idx(j)).events{m}(bidx,3));
                        tf = sum(analysis_in(idx(j)).events{m}(:,3));
                        analysis_in(idx(j)).fractionbound(m) = bf/tf;
                        
                        % remove a value from dwells
                        removeIdx = find(analysis_in(idx(j)).dwells{1}==dt); 
                        analysis_in(idx(j)).dwells{1}(removeIdx(1)) = [];
                    end
                end
            end
        end
        a.events{j} = analysis_in(idx(j)).events;
        a.nevents = a.nevents + sum(analysis_in(idx(j)).nevents);
        a.dwells{1} = [a.dwells{1}; analysis_in(idx(j)).dwells{1}];
        a.dwells{2} = [a.dwells{2}; analysis_in(idx(j)).dwells{2}];
        
        % store dwells CDF
        for l = 1:2
            [cdf1,cdf2,cdf3] = countCDF( a.dwells{l});
            a.dwellsCDF{l} = [cdf1(:),cdf2(:),cdf3(:)];
        end
        
        a.timeToFirst = [a.timeToFirst; analysis_in(idx(j)).timeToFirst];
        a.fractionbound = [a.fractionbound; analysis_in(idx(j)).fractionbound];
    end
    

    if useFitLimits
        dwellLimits = [min(cell2mat(a.frameRate_s)), max(cell2mat(a.nFrames))];
        dwellLimits(2) = dwellLimits(2) * dwellLimits(1); 
    end
    
    if exist('tauGuess', 'var') && ~isempty(tauGuess)
        a.dwellsmle{1} = fitDwells(a.dwells{1}, 2, bootstrap, dwellLimits, tauGuess{1});
        a.dwellsmle{2} = fitDwells(a.dwells{2}, 2, bootstrap, dwellLimits, tauGuess{2});
        a.timeToFirstMLE = fitDwells(a.timeToFirst, 2, bootstrap, dwellLimits, tauGuess{1});
    else
        a.dwellsmle{1} = fitDwells(a.dwells{1}, 2, bootstrap, dwellLimits);
        a.dwellsmle{2} = fitDwells(a.dwells{2}, 2, bootstrap, dwellLimits);
        a.timeToFirstMLE = fitDwells(a.timeToFirst, 2, bootstrap, dwellLimits);
    end
    
    a.unboundmle = [a.dwellsmle{1}.monoExpTau, a.dwellsmle{1}.biExpAmp(1), a.dwellsmle{1}.biExpTau(1), a.dwellsmle{1}.biExpTau(2)];
    a.unboundmle_se = [a.dwellsmle{1}.monoExpSE, a.dwellsmle{1}.biExpAmpSE(1), a.dwellsmle{1}.biExpTauSE(1), a.dwellsmle{1}.biExpTauSE(2)];
    a.boundmle = [a.dwellsmle{2}.monoExpTau, a.dwellsmle{2}.biExpAmp(1), a.dwellsmle{2}.biExpTau(1), a.dwellsmle{2}.biExpTau(2)];
    a.boundmle_se = [a.dwellsmle{2}.monoExpSE, a.dwellsmle{2}.biExpAmpSE(1), a.dwellsmle{2}.biExpTauSE(1), a.dwellsmle{2}.biExpTauSE(2)];
    a.unbound_k = a.dwellsmle{1}.LLR(end);
    a.bound_k = a.dwellsmle{2}.LLR(end);
    
    a.timetofirstmle = [a.timeToFirstMLE.monoExpTau, a.timeToFirstMLE.biExpAmp(1), a.timeToFirstMLE.biExpTau(1), a.timeToFirstMLE.biExpTau(2)];
    a.timetofirstmle_se = [a.timeToFirstMLE.monoExpSE, a.timeToFirstMLE.biExpAmpSE(1), a.timeToFirstMLE.biExpTauSE(1), a.timeToFirstMLE.biExpTauSE(2)];

    [a.fractionbound_fit(1), a.fractionbound_fit(2)] = normfit(a.fractionbound);
    % convert to standard error 
    a.fractionbound_fit(3) = a.fractionbound_fit(2) / sqrt(numel(a.fractionbound)); 
    
     a.field_value = unique_field_values(i);
     a.field_name = field_name;
    
    % store
    if i == 1
        analysis_out = a;
    else
        analysis_out = [analysis_out; a];
    end
end

end