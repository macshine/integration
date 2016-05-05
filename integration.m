function [ci,q,p,z,hc] = integration(data)
%INTEGRATION       Creates time-resolved network topological measures
%
%  [ci,q,p,z,hc] = integration(data);
%
%  This code takes a time-resolved connectivity matrix and 
%  estimates community structure, modularity and the cartographic profile
%  for each region within a region x region x time connectivity matrix.
%  See https://arxiv.org/abs/1511.02976 for more details.
%
%  Requirements: https://sites.google.com/site/bctnet/
%
%  Input:      data     time-series organized in 'nodes x nodes x time' matrix
%
%  Output:     ci       time-resolved community assignment
%              q        time-resolved modularity
%              p        time-resolved participation coefficient
%              z        time-resolved module-degree z-score
%              hc       cartographic profile


    %define variables
    
    [nodes,~,time] = size(data);

    ci = zeros(nodes,time); q = zeros(time,1); p = zeros(nodes,time); z = zeros(nodes,time);

    for t = 1:time
        [ci(:,t),q(t,1)] = modularity_louvain_und_sign(data(:,:,t));
        p(:,t) = participation_coef_sign(data(:,:,t),ci(:,t));
        z(:,t) = module_degree_zscore(data(:,:,t),ci(:,t));
    end

    %cartographic profile

    xbins = [0:0.01:1]; ybins = [5:-.1:-5];
    hc = zeros(size(xbins,2),size(ybins,2),time);
    xNumBins = numel(xbins); yNumBins = numel(ybins);

    for t = 1:time
        Xi = round(interp1(xbins,1:xNumBins,p(:,t),'linear','extrap'));
        Yi = round(interp1(ybins,1:yNumBins,z(:,t),'linear','extrap'));
        Xi = max(min(Xi,xNumBins),1);
        Yi = max(min(Yi,yNumBins),1);
        hc(:,:,t) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
    end

end


