


%%predefine variables
nNodes = 120;
nProbs = 101;
P_out = 0:.01:1;
P_in = 0:.01:1;
mat_id = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4];

ci = zeros(nNodes,nProbs,nProbs);
q = zeros(nProbs,nProbs);
bt = zeros(nNodes,nProbs,nProbs);
wt = zeros(nNodes,nProbs,nProbs);

%%iterate through a range of probabilities for within and between module
%%connectivity

for x = 1:nProbs 
	for y = 1:nProbs

        %%create a graph with fixed a probability of within and between
        %%module edges
        
	mat_out = zeros(nNodes);
	
        for j = 1:nNodes
            mat_out(j,:) = rand(1,nNodes)< P_out(1,y);
            temp = mat_id==mat_id(j,1);
            mat_out(j,temp==1) = rand(1,30)< P_in(1,x);
        end
        
	mat_out = double(mat_out);
	
        %%calculate modular partition, within and between-module connectivity
        %%requires brain connectivity toolbox (https://sites.google.com/site/bctnet/)
        
        [ci(:,x,y),q(x,y)] = community_louvain(mat_out,1); %modular partition
        wt(:,x,y) = module_degree_zscore(mat_out,ci(:,x,y)); %within-module connectivity
        bt(:,x,y) = participation_coef(mat_out,ci(:,x,y)); %between-module connectivity
        
	end
end 

