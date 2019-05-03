%% code for running a simple simulation of within- and between-module connectivity
%
%
%  Plan: create a 4-module network and then populate a proportion of connections
%  either within (P_in) or between (P_out) the modules. Then, calculate a
%  range of different graph theoretical measures to determine how each
%  relates to different topological patterns.

%% probabilities from 0 - 1
P_out = 0:.01:1;
P_in = 0:.01:1;

%% 4 modules, each with 30 nodes
mat_id = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4];

%% predefine variables
nNodes = size(mat_id,1);
mat_out = zeros(nNodes);
nProbs = size(P_out,2);
cig = zeros(nNodes,nProbs,nProbs);
qg = zeros(nProbs,nProbs);
deg = zeros(nNodes,nProbs,nProbs);
bg = zeros(nNodes,nProbs,nProbs);
wg = zeros(nNodes,nProbs,nProbs);
cmg = zeros(nNodes,nNodes,nProbs,nProbs);
ag = zeros(nProbs,nProbs);
gg = zeros(nProbs,nProbs);
lg = zeros(nNodes,nProbs,nProbs);
mat_grp = zeros(nNodes,nNodes,nProbs,nProbs);
nMod = size(unique(mat_id),1);
nSize = nNodes/nMod;


%% loop through probabilities and run graph theoretical analyses

for x = 1:nProbs  
	for y = 1:nProbs

        % for each probability combination (i.e., out and in), fill in that
        % proportion of connections within and between each of the modules
        for j = 1:nNodes
            mat_out(j,:) = rand(1,nNodes)< P_out(1,y);
            temp = mat_id==mat_id(j,1);
            mat_out(j,temp==1) = rand(1,nSize)< P_in(1,x);
        end

        % save each matrix in a larger file for later
        mat_grp(:,:,x,y) = mat_out;
        
        %% requires brain connectivity toolbox (https://sites.google.com/site/bctnet/)
        
        % nodal degree -- the number of connections that each node has
        deg(:,x,y) = degrees_und(mat_out);
        
        % modularity -- estimating the presence of communities
        [cig(:,x,y),qg(x,y)] = community_louvain(mat_out,1);
        
        % quantify the similarity between cig and mat_id using mutual information
        [~,cimi(x,y)] = partition_distance(cig(:,x,y),mat_id);
        
        % cartography -- between ('bg') and within ('wg') module connectivity
        bg(:,x,y) = participation_coef(mat_out,cig(:,x,y));
        wg(:,x,y) = module_degree_zscore(mat_out,cig(:,x,y));
        
        % communicability -- the number of options connecting point A to B
        cmg(:,:,x,y) = expm(mat_out);
        
        % assortavity -- does like connect with like?
        ag(x,y) = assortativity_bin(mat_out,0);
        
        % efficiency -- inverse of characteristic path length
        gg(x,y) = efficiency_bin(mat_out);
        lg(:,x,y) = efficiency_bin(mat_out,1);
        
    end  
    sprintf('%d',x) % gives you a counter so you know how long is left in the simulation
end 




%% make pretty figures
xticklabels = 0:.1:1;
xticks = linspace(1, 101, numel(xticklabels)); 
yticklabels = 0:0.1:1;
yticks = linspace(1, 101, numel(yticklabels));


% degree
figure(1)
imagesc(squeeze(mean(deg,1)))
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Degree')
colorbar
set(gcf,'color','w')

% modularity
figure(2)
imagesc(qg)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Modularity')
colorbar
set(gcf,'color','w')

% similarity between cig and mat_id
figure(3)
imagesc(cimi)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Module recovery')
colorbar
set(gcf,'color','w')

% participation
figure(4)
imagesc(squeeze(mean(bg,1)))
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Participation Coefficient')
colorbar
set(gcf,'color','w')

% module degree Z-score
figure(5)
imagesc(squeeze(mean(wg,1)))
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Module Degree Z-score')
colorbar
set(gcf,'color','w')

% communicability
figure(6)
imagesc(log10(squeeze(mean(mean(cmg,1),2)))) % otherwise the data is too skewed to see
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Communicability')
colorbar
set(gcf,'color','w')

% assortativity
figure(7)
imagesc(squeeze(mean(ag,1)))
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Assortativity')
colorbar
set(gcf,'color','w')

% global efficiency
figure(8)
imagesc(gg)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('Probability of Within-Module Connection')
xlabel('Probability of Between-Module Connection')
title('Global Efficiency')
colorbar
set(gcf,'color','w')
