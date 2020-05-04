% Implementation of different methods for permutation-based knee 
% pathology-pain mapping (KPPM)
% - Lesion and pain data can be simulated or provided.
% - Provided lesions maps have to be spatially aligned before.
% - You can chose between a general linear model (binary logistic regression)
% or Mann-Whitney statistical tests.
% 
% 
% Feel free to use the code under the attached license. Please cite:
% Arthofer et al. An anatomical atlas of the knee and voxel-based 
% knee pathology-pain mapping'
% (C) Copyright 2020, Christoph Arthofer

close all
clear all

% Initialise the random number generator.
rng(42);

% Simulate data or use your own data
create_fake_lesions = false;

% 1... GLM (binary logistic regression)
% 2... Mann-Whitney test (not used at the moment)
stats_test = 1;

% Number of permutations for permutation-based thresholding (currently only
% available with binary logistic regression)
n_perms = 500;

% Output folder for resulting images
output_folder = '~/output/';

if create_fake_lesions
    % Simulate two lesions, one linked to pain, one not linked to pain
    n_imgs = 500;           % number of total images
    max_p = 10;             % pain scale maximum e.g. 0-10
    center1 = [50,50];      % centre of first simulated lesion
    center2 = [120,120];    % centre of second simulated lesion
    patch_dim = 10;      % radius of lesions
    img_size = [160,160];   % size of simulated image

    pain = randi([0,max_p],n_imgs,1);
    mat = zeros([img_size,n_imgs]);
    
    u = sort(unique(pain));
    pain_sub_idx = [];
    for i=1:length(u)
        pain_sub_idx = [pain_sub_idx; randsample(find(pain==i),round(u(i)/max_p*sum(pain==i)))];
    end
    mat(center1(1)-patch_dim:center1(1)+patch_dim, center1(2)-patch_dim:center1(2)+patch_dim,pain_sub_idx) = 1;

    nopain_sub_idx = randperm(n_imgs,numel(pain_sub_idx));
    mat(center2(1)-patch_dim:center2(1)+patch_dim, center2(2)-patch_dim:center2(2)+patch_dim,nopain_sub_idx) = 1;
    
    
    % Write images and design file
    fname = fullfile(output_folder, 'simulated_lesion_pain_design.txt');
    fid = fopen(fname,'w');
    fprintf(fid,'Filename\tpain\tLesionSize\n');
    for i=1:n_imgs
        sim_img = mat2gray(mat(:,:,i));
        imwrite(sim_img,fullfile(output_folder,sprintf('img_%d.png',i)));
        fprintf(fid,'%s\t%d\t%d\n',fname,pain(i),sum(sim_img(:)));
    end
    fclose(fid);

    ref_img = zeros(img_size); % fake template provides a zero backround
    template_name = fullfile(output_folder,'fake_template.png');
    imwrite(mat2gray(ref_img),template_name);

    figure; histogram(pain)
    title('Distribution of pain scores')
    xlabel('Pain score')
    ylabel('# subjects')
    saveas(gcf,fullfile(output_folder, 'distribution_total.png'));
    
    figure; histogram(pain(pain_sub_idx))
    title('Distribution of pain scores in region H');
    xlabel('Pain score')
    ylabel('# subjects')
    saveas(gcf,fullfile(output_folder, 'distribution_region_H.png'));

    figure; histogram(pain(nopain_sub_idx))
    title('Distribution of pain scores in region L');
    xlabel('Pain score')
    ylabel('# subjects')
    saveas(gcf,fullfile(output_folder,'distribution_region_L.png'));

    % should be statistically significant
    idx = zeros(size(pain));
    idx(pain_sub_idx) = 1;
    [h,p] = ttest2(pain(logical(idx)),pain(~logical(idx)))

    % should not be statistically significant
    idx = zeros(size(pain));
    idx(nopain_sub_idx) = 1;
    [h,p] = ttest2(pain(logical(idx)),pain(~logical(idx)))

    mat = reshape(mat,prod(img_size),n_imgs);
    avg = mean(mat,2);
    avg_img = reshape(avg,img_size);
    figure; imagesc(avg_img)
    title('Average image');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    colorbar
    saveas(gcf,fullfile(output_folder,'average.png'));
    
    [tvals, pvals, pvals_perm] = performStats(mat, pain, stats_test, n_perms);
    
    tvals_img = reshape(tvals,img_size);
    figure; imagesc(tvals_img)
    title('t-statistics');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    colorbar
    saveas(gcf,fullfile(output_folder,'t_stats.png'));
    
    pvals_img = reshape(pvals,img_size);
    figure; imagesc(pvals_img)
    title('p-values');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    colorbar
    saveas(gcf,fullfile(output_folder,'p_values.png'));
    
    z = @(p) -sqrt(2) * erfcinv(p*2);
    zvals_img = reshape(abs(z(pvals_img(:))),img_size);
    figure; imagesc(zvals_img)
    title('z-score');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    colorbar
    saveas(gcf,fullfile(output_folder,'z_scores.png'));
else
    % Reference image to get coordinates of reference space; can be a group template (atlas) or random image in reference space
    ref_info = niftiinfo('./anat_knee_template_affine.nii.gz');
    ref_img = niftiread(ref_info);

    % Read data by selecting design file manually
    % Design xlsx file requires the following columns: Filename, pain
    [designFilename,designInputfolder] = uigetfile('~/*.xlsx','Select design file','MultiSelect','off');

    T = readtable(fullfile(designInputfolder,designFilename));
    
    info = niftiinfo(T.Filename{1});
    img_size = info.ImageSize;
    mat = zeros(img_size);
    for s = 1:numel(T.Filename)
        info = niftiinfo(T.Filename{s});
        img = niftiread(info);
        mat(:,:,:,s) = img;
    end
    mat = reshape(mat,prod(img_size),numel(T.Filename));
    
    pain = T.pain(:);
    
    avg = mean(mat,2);
    avg_img = single(reshape(avg,img_size));
    fname = fullfile(output_folder,'average');
    niftiwrite(avg_img,fname,ref_info,'Compressed',true);

    [tvals, pvals, pvals_perm] = performStats(mat, pain, stats_test, n_perms);
    
    tvals_img = single(reshape(tvals,img_size));
    fname = fullfile(output_folder,'t_stats');
	writeImage(tvals_img,fname,ref_info,'Compressed',true);

    pvals_img = single(reshape(pvals,img_size));
    fname = fullfile(output_folder,'p_values');
	writeImage(pvals_img,fname,ref_info);
    
    p_vals_perm_img = single(reshape(pvals_perm,img_size));
    fname = fullfile(output_folder,'p_values_permutation');
    writeImage( 1-p_vals_perm_img,fname,ref_info,'Compressed',true);

    z = @(p) -sqrt(2) * erfcinv(p*2);
    zvals_img = single(reshape(abs(z(pvals_img(:))),img_size));
    fname = fullfile(output_folder,'z_scores');
	writeImage(zvals_img,fname,ref_info,'Compressed',true);
end

warning('on','all')



function [tvals, pvals, pvals_perm] = performStats(mat, pain, stats_test, nperms)
% Function performing the permutation-based statistical tests.
% Input:
%   mat:            N x P matrix containing spatially aligned images; N voxels x P images
%   pain:           column vector of pain scores
%   stats_test:     1...binary logistic regression (with permutations), 2...Mann-Whitney (without permutations)
%   nperms:         number of permutations
% Output:
%   tvals:          column vector of t-statistic; N x 1
%   pvals:          column vector of p-values; N x 1
%   pvals_perm:     column vector of p-values after permutations; N x 1

    avg = mean(mat,2);
    size_avg = size(avg(:));
    size_mat = size(mat);
    
    pvals = ones(size_avg);
    
    ind = find(avg);
    nind = length(ind);
    sprintf('Nind: %d',nind);
    
    if stats_test == 1
        % Binary Logistic Regression
        sprintf('Binary logistic regression')
        
        tvals = zeros(size_avg);
        pvals_perm = zeros(size_avg);
    
        X = pain(:);
        tvals_temp = zeros(nind,1);
        pvals_temp = zeros(nind,1);
        parfor i = 1:nind
            y = mat(ind(i),:)';
            mdl = fitglm(X,y,'Distribution','binomial');
            tvals_temp(i) = mdl.Coefficients.tStat(2);
            pvals_temp(i) = mdl.Coefficients.pValue(2);
        end
        tvals(ind) = tvals_temp;
        pvals(ind) = pvals_temp;
        
        if nperms > 1
        sprintf('Number of permutations: %d',nperms)
            for j = 1:nperms
                sprintf('Permutation: %d',j)
                X = pain(:);
                rand_idx = randperm(size_mat(2));
                t_vals_perm_temp = zeros(nind,1);
                t_vals_perm = zeros(size_avg);
                parfor i = 1:nind
                    y = mat(ind(i),:)';
                    mdl = fitglm(X,y(rand_idx),'Distribution','binomial');
                    t_vals_perm_temp(i) = mdl.Coefficients.tStat(2);
                end
                t_vals_perm(ind) = t_vals_perm_temp;
                idxs = tvals <= t_vals_perm;
                pvals_perm(idxs) = pvals_perm(idxs) + 1;
            end
            pvals_perm = pvals_perm ./ nperms;
        end
    elseif stats_test == 2
        sprintf('Currently not in use!');
%         % Mann-Whitney test
%         sprintf('Mann-Whitney test.... \nOnly a p-value map will be provided...')
%         for i = 1:nind
%             lesion_idx = mat(ind(i),:) > 0;
%             x = pain(lesion_idx);
%             y = pain(~lesion_idx);
%             [p,h,stats] = ranksum(x,y);
%             pvals(ind(i)) = p;
%         end
    end
end









