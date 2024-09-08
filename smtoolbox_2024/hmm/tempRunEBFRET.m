[file,path] = uigetfile('.mat');
load([path,file]);
disp(file)

%% Convert each selected trace to a cell 
if length(data.time{1}(:,end)) >= length(data.time{2}(:,end))
    disp('crop')
    remove_idx = removeSameOffEvent(data.rois);
    rois = data.rois(:,2);
    rois(remove_idx,:) = [];
else
    rois = data.rois(:,2);
end
rois = rois(vertcat(rois.status)>0);
nrois = length(rois);

%% 
traces = cell(nrois,1); 
for i = 1:nrois
    traces{i} = rois(i).timeSeries;
end

%%
% u = guess_prior(x, K, spread, strength)
u0 = guess_prior(traces, 2); 

%%
[u, L, vb, vit, phi] = ebayes(traces, u0)