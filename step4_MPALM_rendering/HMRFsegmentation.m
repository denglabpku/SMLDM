function [HMRFseg, settings] = HMRFsegmentation(raw_table, xylim, pre_nucleus_mask, nclust, dxy, sigma_render, beta0,ImPALM)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明

    % ----------------------- parameters settings -------------------------
    maxiter = 200;
    mineps = 10^(-7);
    varfixed = true;
    inforce = true;
    if nargin < 4
        dxy = 20;
        sigma_render = 50;
        nclust = 7;
        beta0 = 0.1;
    end

    % ------------------------ nucleus mask ---------------------------
    % Xmin = xylim(1); Xmax = xylim(2); Ymin = xylim(3); Ymax = xylim(4);
    % delt_X=Xmax-Xmin;delt_Y=Ymax-Ymin;
    % delta=max(delt_X,delt_Y);
    % num_pix = round(delta/dxy);     
    % ImPALM = uint16(Render2(raw_table.Xpos, raw_table.Ypos, xylim, dxy, sigma_render, true, false));
   
    if sum(pre_nucleus_mask(:)) == 0
        bw0 = stdfilt(ImPALM, true(5));
        bw0 = imgaussfilt(bw0, 3);
        bw = imbinarize(bw0, 15);
        bw = imfill(bw, 'holes');
        ccomp = bwconncomp(bw);
        nucleus_mask = zeros(size(ImPALM));
        [~, i] = max(cellfun(@length, ccomp.PixelIdxList));
        nucleus_mask(ccomp.PixelIdxList{1, i}) = 1;
        nucleus_mask = logical(nucleus_mask');
    else
        nucleus_mask = pre_nucleus_mask;
    end
    % imagesc(nucleus_mask);
    
    % ---------- Hidden Markov Random Field (HMRF) model --------------
 
    X = size(nucleus_mask, 1); Y = size(nucleus_mask, 2);
    mask = nucleus_mask(:);

    % img = Render2(raw_table.Xpos, raw_table.Ypos, xylim, num_pix, sigma_render, true, false);
    img = ImPALM;
    
    % initiation
    img = img'; 
    img_class = zeros(size(mask)); % from 1 to nclust
    normf = quantile(img(mask), 0.95); % norm factor not contribute significant difference
    img = img(:)/normf; 
    qu = quantile(img(mask), (1:(nclust-1))/nclust);
    for i = 1:(nclust-1)
        img_class((img>qu(i))&mask) = i;
    end
    img_class(mask) = img_class(mask) + 1;
    beta = diag(ones(1, nclust)*beta0);
    
    mu = zeros(nclust, 1);
    sigma = ones(nclust, 1);
    
    for i = 1:nclust
        mu(i) = mean(img(img_class==i));
        sigma(i) = std(img(img_class==i));
    end
    if varfixed
        sigma = ones(nclust, 1);
        sigma = std(img(mask)-mu(img_class(mask)))*sigma;
    end
    loglik = zeros(nclust, 1);
    counter = 0;
    criterium = true;
    
    % interation
    while criterium
        counter = counter + 1;
        disp(['Iteration: ', num2str(counter)]);
    
        for xidx = 1:X
            for yidx = 1:Y
                idx = (yidx-1)*X + xidx;
                if mask(idx) 
                    for i = 1:nclust
                        loglik(i) = -0.5*(img(idx)-mu(i))^2/sigma(i)/sigma(i);
                    end
                    if ~(xidx == 1)
                        nid = (yidx-1)*X+xidx-1;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    if ~(xidx == X)
                        nid = (yidx-1)*X+xidx+1;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    if ~(yidx == 1)
                        nid = (yidx-2)*X+xidx;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    if ~(yidx == Y)
                        nid = yidx*X+xidx;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    [~, temp] = max(loglik);
                    img_class(idx) = temp;
                end
            end
        end
    
        temp_nclust = nclust;
        for i = nclust:-1:1
            if sum(img_class == i) == 0
                disp(["class ", num2str(i), " removed."]);
                img_class(img_class>i) = img_class(img_class>i)-1;
                temp_nclust = temp_nclust - 1;
            end
        end
        oldmu = mu(1:temp_nclust);
        mu = zeros(temp_nclust, 1);
        sigma = ones(temp_nclust, 1);
        for i = 1:temp_nclust
            mu(i) = mean(img(img_class==i));
            sigma(i) = std(img(img_class==i));
        end
        if varfixed
            sigma = ones(nclust, 1);
            sigma = std(img(mask)-mu(img_class(mask)))*sigma;
        end
        [~,~,ic] = unique(img_class);
        a_counts = accumarray(ic,1);
        sigma(sigma==0) = 1e-6;
        disp("Class number: ");
        disp(num2str(a_counts));
        disp("Class mean:");
        disp(num2str(mu));
        disp("Class sigma:");
        disp(num2str(sigma));
    
        if (counter == maxiter)||(sum((mu-oldmu).^2)<mineps)
            criterium = false;
        end
    
        if inforce
            while temp_nclust < nclust
                criterium = true;
                disp(["inforce nclust: ", num2str(nclust)]);
                if ~varfixed
                    [~, w] = max(sigma);
                else
                    [~,~,ic] = unique(img_class);
                    a_counts = accumarray(ic,1);
                    [~, w] = max(a_counts(2:end));
                end
                img_class(img_class>w) = img_class(img_class>w) + 1;
                nn = sum(img_class == w);
                disp("Class with max sigma has pixel number: ", num2str(nn));
                img_class((img_class==w)&(img>mu(w))) = img_class((img_class==w)&(img>mu(w))) + 1;
                temp_nclust = temp_nclust + 1;
                mu = zeros(temp_nclust, 1);
                sigma = ones(temp_nclust, 1);
                for i = 1:temp_nclust
                    mu(i) = mean(img(img_class==i));
                    sigma(i) = std(img(img_class==i));
                end
                if varfixed
                    sigma = ones(nclust, 1);
                    sigma = std(img(mask)-mu(img_class(mask)))*sigma;
                end
                sigma(sigma==0) = 1e-6;
            end
        end
    end
    disp(['Stop at iteration: ', num2str(counter)]);
    
    % --- ZW modified to keep all initial classification --------------- %
    % if sum(pre_nucleus_mask(:)) == 0
    %     temp_region = img_class > 1;
    %     temp_region = reshape(temp_region, size(nucleus_mask));
    %     nucleus_mask = imfill(temp_region, "holes");
    %     ccomp = bwconncomp(nucleus_mask>0);
    %     nucleus_mask = false(size(nucleus_mask));
    %     [~, i] = max(cellfun(@length, ccomp.PixelIdxList));
    %     nucleus_mask(ccomp.PixelIdxList{1, i}) = 1;
    % end
    % mask = nucleus_mask(:);
    % img_class(~mask) = 0;
    % [~,~,ic] = unique(img_class);
    % a_counts = accumarray(ic,1); a_counts = a_counts(2:end);

    % save results and parameters
    HMRFseg.img = img;
    HMRFseg.img_class = img_class;
    HMRFseg.nucleus_mask = nucleus_mask; % HMRFseg.mask = mask; 
    HMRFseg.mu = mu;
    HMRFseg.sigma = sigma;
    HMRFseg.a_counts = a_counts;

    settings.dxy = dxy; settings.maxiter = maxiter; settings.normf = normf;
    settings.mineps = mineps; settings.nclust = nclust; settings.beta0 = beta0;
    settings.varfixed = varfixed; settings.inforce = inforce; settings.sigma_render = sigma_render;
end