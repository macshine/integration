%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Step 1: Set-Up

%addpath to scripts
addpath('/path/to/file/2014_04_05 BCT/') % https://sites.google.com/site/bctnet/

%load preprocessed data (rows = regions & columns = time)
data = dlmread('data.tsv','\t',1,1);

%load ID matrices - Gordon et al. network id (http://www.ncbi.nlm.nih.gov/pubmed/25316338)
load('/path/to/file/id1.mat');

%identify variable sizes
[nNodes,nTime] = size(data);


%%Step 2: Functional Connectivity

%time-averaged codnnectivity matrix
stat_avg = corr(data');

%time-resolved connectivity - Multiplication of Temporal Derivatives (MTD)
td = diff(data');
data_std = std(td);

for n = 1:nNodes
  td(:,n) = td(:,n) / data_std(1,n);
end

raw_fc = bsxfun(@times,permute(td,[1,3,2]),permute(td,[1,2,3]));

%Simple moving average of MTD
w = 14; % window length = 14 TRs (~10 seconds using 0.72s TR data)
sma_filter = 1/w*ones(w,1);
sma = zeros(nTime,nNodes,nNodes);

for j = 1:nNodes
  for k = 1:nNodes
    sma(2:end,j,k) = filter(sma_filter,1,raw_fc(:,j,k));
  end
end

sma = permute(sma,[2,3,1]);

%time-averaged connectivity matrix
dyn_avg = nanmean(sma,3);
dyn_z = weight_conversion(dyn_avg,'normalize'); %normalize


%% Step 3: Graph Theoretical Measures

%Modularity
ci = zeros(nNodes,nTime);
q = zeros(nTime,1);

for t = 1:nTime
  [ci(:,t),q(t,1)] = modularity_louvain_und_sign(sma(:,:,t));
end

q_avg = nanmean(q);


%Module Degree Z-score
mod_deg_z = zeros(nNodes,nTime);

for t = 1:nTime
  mod_deg_z(:,t) = module_degree_zscore(sma(:,:,t),ci(:,t),0);
end

Z_avg = nanmean(mod_deg_z,2);
Z_std = nanstd(mod_deg_z(:));


%Participation index
P = zeros(nNodes,nTime);

for t = 1:nTime
  P(:,t) = participation_coef_sign(sma(:,:,t),ci(:,t));
end

P_avg = nanmean(P,2);
P_std = nanstd(P(:));


% Step 4: Cartographic Analysis (http://www.nature.com/nature/journal/v433/n7028/full/nature03288.html)


% Step 5: 2-dimensional histogram cloud

xbins = [0:0.01:1.0]; ybins = [5:-.1:-5]; % 100 x 100 2d histogram
  
hist_cloud = zeros(size(xbins,2),size(ybins,2),nTime); %predefine for speed

xNumBins = numel(xbins); yNumBins = numel(ybins);

for t = 1:nTime
  Xi = round(interp1(xbins, 1:xNumBins, P(:,t), 'linear', 'extrap') );
  Yi = round(interp1(ybins, 1:yNumBins, mod_deg_z(:,t), 'linear', 'extrap') );
  Xi = max( min(Xi,xNumBins), 1);
  Yi = max( min(Yi,yNumBins), 1);
  hist_cloud(:,:,t) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
end

% Step 6: K-means analysis

idx = kmeans(reshape(hist_cloud,xNumBins * yNumBins,nTime),2);


