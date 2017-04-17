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

%calculate temporal derivative
td = diff(data);

%functional coupling score
fc = bsxfun(@times,permute(td,[1,3,2]),permute(td,[1,2,3]));

%simple moving average (need to define window length)

window = 15;

mtd_filter = 1/window*ones(window,1);
mtd = zeros(nTime,nNodes,nNodes);

for j = 1:nNodes
    for k = 1:nNodes
        mtd(2:end,j,k) = filter(mtd_filter,1,fc(:,j,k));
    end
end

mtd(1:round(window/2),:,:) = [];
mtd(round(nTime-window):nTime,:,:) = 0;
mtd = permute(mtd,[2,3,1]);
mtd(:,:,1) = mtd(:,:,2);

%time-averaged connectivity matrix
dyn_avg = nanmean(mtd,3);
dyn_z = weight_conversion(dyn_avg,'normalize'); %normalize



%% Step 3: Graph Theoretical Measures

%Modularity
ci = zeros(nNodes,nTime);
q = zeros(nTime,1);

for t = 1:nTime
  [ci(:,t),q(t,1)] = modularity_louvain_und_sign(sma(:,:,t)); %%this step should be run multiple times in order to obtain consensus
end


%Module Degree Z-score (WT)
WT = zeros(nNodes,nTime);

for t = 1:nTime
  WT(:,t) = module_degree_zscore(sma(:,:,t),ci(:,t),0);
end

%Participation index (BT)
BT = zeros(nNodes,nTime);

for t = 1:nTime
  BT(:,t) = participation_coef_sign(sma(:,:,t),ci(:,t));
end



% Step 4: 2-dimensional Cartographic Profile (CP)

xbins = [0:0.01:1.0]; ybins = [5:-.1:-5]; % 100 x 100 2d histogram
CP = zeros(size(xbins,2),size(ybins,2),nTime);
xNumBins = numel(xbins); yNumBins = numel(ybins);

for t = 1:nTime
  Xi = round(interp1(xbins, 1:xNumBins, BT(:,t), 'linear', 'extrap') );
  Yi = round(interp1(ybins, 1:yNumBins, WT(:,t), 'linear', 'extrap') );
  Xi = max( min(Xi,xNumBins), 1);
  Yi = max( min(Yi,yNumBins), 1);
  CP(:,:,t) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
end


% Step 5: K-means analysis

idx = kmeans(reshape(hist_cloud,xNumBins * yNumBins,nTime)',2); %%also worth increasing k to determine stability of clustering


