function [ci_new,F] = hungarian1(ci,beta)

    %unique modules per timepoint
    [nNodes,nTime] = size(ci);
    nDer = nTime - 1;
    
%     dyn_mod = zeros(nMod,nTime);

    for t = 1:nTime
        temp = tabulate(ci(:,t));
        dyn_mod(1:size(temp,1),t) = temp(:,1);
    end    
    
    %number of modules per window
    number_mod = zeros(nTime,1);

    for t = 1:nTime
        number_mod(t,1) = nnz(dyn_mod(:,t));
    end

    ignore = double(number_mod>5);
    nMod = max(max(ci(:,ignore==0)));
    dyn_mod(nMod+1:end,:) = [];
    
    %1-of-k encoding
    encode = zeros(nNodes,nMod,nTime);

    for t = 1:nTime
        if ignore(t,1)==0
            C = ci(:,t);
            R = 1:numel(C);
            A = zeros(numel(C),max(C));
            A(sub2ind(size(A),R',C)) = 1;
            encode(:,1:size(A,2),t) = A;
        else
            encode(:,:,t) = 0;
        end
    end

    sprintf('%s','1 of k encoding')
    

    %dice between consecutive time points
        %this is usually bimodal
    dice_coef_encode = zeros(nMod,nMod,nTime);

    for t = 1:nDer
        dice_coef_encode(:,:,t) = bsxfun(@corr,encode(:,:,t),encode(:,:,t+1));
    end

    sprintf('%s','similarity')
    

    %threshold & cost
    cost = 1/(double(dice_coef_encode>beta));


    %hungarian algorithm
    assignment = zeros(nTime,nMod);

    for t = 1:nTime
        [assignment(t,:),~] = munkres(cost(:,:,t));
    end

    sprintf('%s','munkres')
    
    
    
    
    %hungarian un-twisting 
    dyn_mod2 = zeros(nMod,nTime);  

    for t = 1:nTime
        for k = 1:nMod
            if number_mod(t,1) < k
                dyn_mod2(k,t) = NaN;
            end
        end
    end

    dyn_mod2(:,1,:) = dyn_mod(:,1,:); %dyn_mod2 starting with dyn_mod's first assignment


    % recoding
    tally = max(dyn_mod2(:,1));            

    for w = 2:nTime-1
        for k = 1:nMod
            for l = 1:nMod
                if dyn_mod2(k,w-1) == 0
                    dyn_mod2(k,w-1) = tally+1;
                    tally = tally+1;
                end

                if assignment(w-1,k) == l
                    dyn_mod2(l,w) = dyn_mod2(k,w-1);
                end
            end
        end
    end

    
%     % number of timepoints that a module exists for
%     count_mod = zeros(nMod,1);
% 
%     for x = 1:max(dyn_mod2(:))
%         count_mod(x,1) = countmember(x,dyn_mod2(:));
%     end
% 
% 
%     % parcel size of each module over time
%     count_mod_size = zeros(nMod,nTime);
% 
%     for t = 1:nTime
%         for y = 1:nMod
%             count_mod_size(y,t) = countmember(y,ci(:,t));
%         end
%     end
% 
%     for t = 1:nTime
%         count_mod_unique(1:size(unique(count_mod_size(:,t))),t) = unique(count_mod_size(:,t));
%     end
% 
%     count_mod_unique(1,:) = [];


    % recode the net_thr matrix using new module names (~mucha)
    

    ci_temp = zeros(nNodes,nTime);

    for t = 1:nTime
        for j = 1:nNodes
            for k = 1:nMod
                if ci(j,t) == dyn_mod(k,t)
                    ci_temp(j,t) = dyn_mod2(k,t);
                end
            end
        end
    end

    sprintf('%s','rename')
    

    %%temporal sorting

    %ci_new_name_1_of_k
    encode2 = zeros(nNodes,tally,nTime);

    for t = 1:nTime
        for j = 1:nNodes
            for p = 1:tally
                if ci_temp(j,t) == p
                    encode2(j,p,t) = 1;
                end
            end
        end
    end


    %ci_signature
    encode2_sum = nanmean(encode2,3);
    encode2_perc = encode2_sum / tally;


    %ci_correlation matrix -- can we use a better similarity metric here?
    corr_sig = corr(encode2_perc);

    
    %use affinity propogation to cluster
    [ap,~,~,~] = apcluster(corr_sig,nanmin(corr_sig(:)));

    sprintf('%s','affinity propagation')
    
    
    %identity of clusters estimated by affinity propogation
    unique_ap = unique(ap);

    
    %number of clusters estimated by affinity propogation
    count_ap = nnz(unique(ap));

    for x = 1:size(unique_ap,1)
        freq_ap(x,1) = sum(ap==unique_ap(x));
    end
    
    
    %rename ci_new_name according to ap clusters
    dyn_mod3 = dyn_mod2;

    for t = 1:nTime
        for k = 1:nMod
            for p = 1:tally
                if dyn_mod2(k,t) == p
                    dyn_mod3(k,t) = ap(p,1);
                end
            end
        end
    end

    
    % rename parcels according to affinity clustering
    ci_temp2 = zeros(nNodes,nTime);
    ci_new = zeros(nNodes,nTime);
    

    for t = 1:nTime
        for j = 1:nNodes
            for k = 1:nMod
                if ci_temp(j,t) == dyn_mod2(k,t)
                    ci_temp2(j,t) = dyn_mod3(k,t);
                end
            end
        end
    end
    
    
    %% rename parcels according to order of appearance
    
    for t = 1:nTime
        for j = 1:nNodes
            for x = 1:size(unique_ap,1)
                y = unique_ap(x);
                if ci_temp2(j,t) == y
                    ci_new(j,t) = find(unique_ap==y);
                end
            end
        end
    end
        
    F = flexibility(ci_new');
    
    sprintf('%s','flexibility')
    
    
end


    