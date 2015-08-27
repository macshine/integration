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
W0 = cov(data'); % this will be used to create stationary (VAR) data

%time-resolved connectivity - Multiplication of Temporal Derivatives (MTD)
td = diff(data');
data_std = std(td);

for n = 1:nNodes
  td(:,n) = td(:,n) / data_std(1,n);
end

raw_fc = bsxfun(@times,permute(td,[1,3,2]),permute(td,[1,2,3]));

%Simple moving average of MTD
w = 10; % window length = 10 TRs
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


%Degeneracy
ci_deg = zeros(nNodes,nTime);

for t = 1:nTime
  ci_deg(:,t) = modularity_finetune_und_sign(sma(:,:,t),'sta',ci(:,t));
end

%Module Degree Z-score
mod_deg_z = zeros(nNodes,nTime);

for t = 1:nTime
  mod_deg_z(:,t) = module_degree_zscore(sma(:,:,t),ci_deg(:,t),0);
end

Z_avg = nanmean(mod_deg_z,2);
Z_std = nanstd(mod_deg_z(:));


%Participation index
P = zeros(nNodes,nTime);

for t = 1:nTime
  P(:,t) = participation_coef_sign(sma(:,:,t),ci_deg(:,t));
end

P_avg = nanmean(P,2);
P_std = nanstd(P(:));


% Step 4: Cartographic Analysis (http://www.nature.com/nature/journal/v433/n7028/full/nature03288.html)

Z_P1 = zeros(nNodes,nTime);
Z_P2 = zeros(nNodes,nTime);
Z_P3 = zeros(nNodes,nTime);
Z_P4 = zeros(nNodes,nTime);
Z_P5 = zeros(nNodes,nTime);
Z_P6 = zeros(nNodes,nTime);
Z_P7 = zeros(nNodes,nTime);

for t = 1:nTime
  for j = 1:nNodes
    if Z_ci(j,t) < 2.5 & P(j,t) < 0.05 
      Z_P1(j,t) = 1; % ultra-peripheral nodes
    elseif Z_ci(j,t) < 2.5 & P(j,t) >= 0.05 & P(j,t) < 0.62 
      Z_P2(j,t) = 1; % peripheral nodes
    elseif Z_ci(j,t) < 2.5 & P(j,t) >= 0.62 & P(j,t) < 0.8 
      Z_P3(j,t) = 1; % connector nodes
    elseif Z_ci(j,t) < 2.5 & P(j,t) >= 0.8
      Z_P4(j,t) = 1; % kinless nodes
    elseif Z_ci(j,t) >= 2.5 & P(j,t) < 0.3
      Z_P5(j,t) = 1; % provincial hubs
    elseif Z_ci(j,t) >= 2.5 & P(j,t) >= 0.3 & P(j,t) < 0.75
      Z_P6(j,t) = 1; % connector hubs
    elseif Z_ci(j,t) >= 2.5 & P(j,t) >= 0.75
      Z_P7(j,t) = 1; % kinless hubs
    end
  end
end

