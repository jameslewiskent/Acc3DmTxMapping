function [] = SimulateUndersamplingofSyntheticBodyImages(Recon_Type,sx,sy,sz,calib,niters_array,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)
load('SyntheticBodyImages.mat','Relative_Rx','Ref_Shims_mTx','Prep_Shims_mTx','image_mask_full');

kernel = [5 5]; % kernel size
eig_thresh = 0.02; % Threshold for picking singular vectors of the calibration matrix (relative to largest singular value)
eig_thresh2 = 0.95; % Threshold of eigenvector decomposition in image space
calibsize = [sx sy]; % ESPIRiT calibration size

location = ['Data',filesep,'Synthetic Body Simulation Results',filesep];
Folder1 = [location,'AcceleratedData',filesep];
Folder2 = [location,'ReconData',filesep];

% Centre body in FOV
image_mask_full = circshift(image_mask_full,15,1);
Relative_Rx = circshift(Relative_Rx,15,1);
Ref_Shims_mTx = circshift(Ref_Shims_mTx,15,1);
Prep_Shims_mTx = circshift(Prep_Shims_mTx,15,1);
%% Crop out arms (optional)
%figure(); imagesc(imtile(image_mask_full))
image_mask_full(:,1:8,:,:)= 0; image_mask_full(:,170:178,:,:)= 0;
%figure(); imagesc(imtile(image_mask_full))
Relative_Rx = Relative_Rx.*image_mask_full;
Ref_Shims_mTx = Ref_Shims_mTx.*image_mask_full;
Prep_Shims_mTx = Prep_Shims_mTx.*image_mask_full;

%% FFT and crop k-space to desired size (sx, sy)
kdata_Relative = fft2c(Relative_Rx);
kdata_Ref_Shims_mTx = fft2c(Ref_Shims_mTx);
kdata_Prep_Shims_mTx = fft2c(Prep_Shims_mTx);

% Crop kspace to sx * sy
min_ind_kx = ceil(round((size(kdata_Relative,1)+1)/2) - sx/2);
min_ind_ky = ceil(round((size(kdata_Relative,2)+1)/2) - sy/2);
k_Relativec = 1e6*kdata_Relative(min_ind_kx:min_ind_kx+sx-1,min_ind_ky:min_ind_ky+sy-1,:,:);
k_Ref_Shims_mTxc = 1e6*kdata_Ref_Shims_mTx(min_ind_kx:min_ind_kx+sx-1,min_ind_ky:min_ind_ky+sy-1,:,:);
k_Prep_Shims_mTxc = 1e6*kdata_Prep_Shims_mTx(min_ind_kx:min_ind_kx+sx-1,min_ind_ky:min_ind_ky+sy-1,:,:);

% define noise levels
std1 = max(abs(k_Relativec),[],'all')./db2mag(pSNR); std2 = max(abs(k_Ref_Shims_mTxc),[],'all')./db2mag(pSNR); std3 = max(abs(k_Prep_Shims_mTxc),[],'all')./db2mag(pSNR);
%% Under-sample images using under-sampling masks, then reconstruct k-space using admm
tic
for iter_n  = 1:size(niters_array,2)
    niters = niters_array(iter_n);
    filename = ['Simulated_',Recon_Type,'Recon_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats)];
    
    % pre-allocate
    recon_rel = zeros([sx,sy,8,size(k_Relativec,4),NRepeats,3,size(accelerations,2)]);
    if strcmp(Recon_Type,'SENSE') ||  strcmp(Recon_Type,'GRAPPA')
        refUS = zeros([sx,sy,8,size(k_Ref_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]);
        prepUS = zeros([sx,sy,8,size(k_Prep_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]);
    else
        recon_ref = zeros([sx,sy,8,size(k_Ref_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]);
        recon_prep = zeros([sx,sy,8,size(k_Prep_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]);
    end
    
    if exist([Folder1,filename,'.mat'],'file') == 2
        disp(['File ', filename,'.mat', ' already exists, loading in data'])
        load([Folder1,filename,'.mat']);
    else
        parfor accel_ind = 1:size(accelerations,2)
            kdata_Relativec_acc = zeros(sx,sy,8,8);
            kdata_Ref_Shims_mTxc_acc = zeros(sx,sy,8,2);
            kdata_Prep_Shims_mTxc_acc = zeros(sx,sy,8,2);
            for same_masks = 2 %1:3 Code used for testing different alternative masking strategies
                for repeat_n = 1:NRepeats
                    % Corrupt with noise
                    k_Relativecn = k_Relativec + std1*randn(size(k_Relativec)) + 1i*std1*randn(size(k_Relativec));
                    k_Ref_Shims_mTxcn = k_Ref_Shims_mTxc + std2*randn(size(k_Ref_Shims_mTxc)) + 1i*std2*randn(size(k_Ref_Shims_mTxc));
                    k_Prep_Shims_mTxcn  = k_Prep_Shims_mTxc  + std3*randn(size(k_Prep_Shims_mTxc)) + 1i*std3*randn(size(k_Prep_Shims_mTxc));
                    
                    % Controls same/different masks for Tx modes and
                    % ref/prepared images
                    if same_masks == 1 % all Tx modes, all ref/prep images use same masks
                        rel_mask_n = repeat_n;
                        ref_mask_n = repeat_n;
                        prep_mask_n = repeat_n;
                        % apply masks
                        kdata_Relativec_acc = bsxfun(@times,k_Relativecn,masks(:,:,accel_ind,rel_mask_n));
                        kdata_Ref_Shims_mTxc_acc = bsxfun(@times,k_Ref_Shims_mTxcn,masks(:,:,accel_ind,ref_mask_n));
                        kdata_Prep_Shims_mTxc_acc = bsxfun(@times,k_Prep_Shims_mTxcn,masks(:,:,accel_ind,prep_mask_n));
                    elseif same_masks == 2 % Tx modes use different masks, ref/prep images use same masks (paired)
                        n_refprep_masks = size(k_Prep_Shims_mTxc,4); n_rel_masks = size(k_Relativec,4); count = 1; total_per_rep = n_refprep_masks + n_rel_masks;
                        rep_base = ((repeat_n-1)*total_per_rep)+1;
                        rel_mask_n = (rep_base:rep_base+n_rel_masks-1);
                        ref_mask_n = (rep_base+n_rel_masks:rep_base+n_rel_masks+n_refprep_masks-1);
                        prep_mask_n = ref_mask_n; % reference image uses same mask as prepared
                        % apply masks
                        kdata_Relativec_acc = bsxfun(@times,k_Relativecn,masks(:,:,accel_ind,rel_mask_n));
                        kdata_Ref_Shims_mTxc_acc = bsxfun(@times,k_Ref_Shims_mTxcn,masks(:,:,accel_ind,ref_mask_n));
                        kdata_Prep_Shims_mTxc_acc = bsxfun(@times,k_Prep_Shims_mTxcn,masks(:,:,accel_ind,prep_mask_n));
                    elseif same_masks == 3 % all Tx modes, all ref/prep images use different masks
                        n_refprep_masks = size(k_Prep_Shims_mTxc,4); n_rel_masks = size(k_Relativec,4); count = 1; total_per_rep = 2*n_refprep_masks + n_rel_masks;
                        rep_base = ((repeat_n-1)*total_per_rep)+1;
                        rel_mask_n = (rep_base:rep_base+n_rel_masks-1);
                        ref_mask_n = (rep_base+n_rel_masks:rep_base+n_rel_masks+n_refprep_masks-1);
                        prep_mask_n = (rep_base+n_rel_masks+n_refprep_masks:rep_base+n_rel_masks+2*n_refprep_masks-1);
                        % apply masks
                        kdata_Relativec_acc = bsxfun(@times,k_Relativecn,masks(:,:,accel_ind,rel_mask_n));
                        kdata_Ref_Shims_mTxc_acc = bsxfun(@times,k_Ref_Shims_mTxcn,masks(:,:,accel_ind,ref_mask_n));
                        kdata_Prep_Shims_mTxc_acc = bsxfun(@times,k_Prep_Shims_mTxcn,masks(:,:,accel_ind,prep_mask_n));
                    end
                    
                    if  strcmp(Recon_Type,'jTxLR')
                        joint_recon = admm_txlr(double(cat(4,kdata_Relativec_acc,kdata_Ref_Shims_mTxc_acc,kdata_Prep_Shims_mTxc_acc)), kernel, niters, [50 50]);
                        recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind) = joint_recon(:,:,:,1:8);
                        recon_ref(:,:,:,:,repeat_n,same_masks,accel_ind) = joint_recon(:,:,:,9:10);
                        recon_prep(:,:,:,:,repeat_n,same_masks,accel_ind) = joint_recon(:,:,:,11:12);
                    else
                        % Joint recon of relative images by concatenating along Tx dimension
                        recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind) = admm_txlr(double(kdata_Relativec_acc), kernel, niters, [50 50]);
                        if  strcmp(Recon_Type,'SENSE') ||  strcmp(Recon_Type,'GRAPPA')
                            refUS(:,:,:,:,repeat_n,same_masks,accel_ind) = kdata_Ref_Shims_mTxc_acc; % reference and prepared images still undersampled
                            prepUS(:,:,:,:,repeat_n,same_masks,accel_ind) = kdata_Prep_Shims_mTxc_acc;
                        else
                            recon_ref(:,:,:,:,repeat_n,same_masks,accel_ind) = admm_txlr(double(kdata_Ref_Shims_mTxc_acc), kernel, niters, [50 50]);
                            recon_prep(:,:,:,:,repeat_n,same_masks,accel_ind) = admm_txlr(double(kdata_Prep_Shims_mTxc_acc), kernel, niters, [50 50]);
                        end
                    end
                end
            end
            disp(['Acceleration Factor of ',num2str(accelerations(accel_ind)), ' completed ', datestr(now, 'dd/mm/yy-HH:MM')])
        end
        if  strcmp(Recon_Type,'SENSE') ||  strcmp(Recon_Type,'GRAPPA')
            save([Folder1,filename,'.mat'],'recon_rel','refUS','prepUS')
        else
            save([Folder1,filename,'.mat'],'recon_rel','recon_ref','recon_prep')
        end
    end
end
toc

% Check nmrse
%nrmse(k_Relativecn,k_Relativec) % Fully sampled noisey vs noiseless
%nrmse(recon_rel(:,:,:,:,1,2,1),k_Relativec) % Fully sampled noisey (been through TxLR) vs noiseless
%nrmse(recon_rel(:,:,:,:,1,2,3),k_Relativec) % Undersampled noisey (been through TxLR) vs noiseless
%nrmse(recon_ref(:,:,:,:,1,2,3),k_Ref_Shims_mTxc) % Undersampled noisey (been through TxLR) vs noiseless
%nrmse(recon_prep(:,:,:,:,1,2,3),k_Prep_Shims_mTxc) % Undersampled noisey (been through TxLR) vs noiseless

clear k_Relativecn k_Ref_Shims_mTxcn k_Prep_Shims_mTxcn
%% Define mask
% reduce size of image mask by FFT and cropping in k-space
k_mask = fft2c(image_mask_full);
min_ind_kz1 = round(round((size(k_mask,1)+1)/2) - sz(1)/2);
min_ind_kz2 = round(round((size(k_mask,2)+1)/2) - sz(2)/2);
k_mask = k_mask(min_ind_kz1:min_ind_kz1+sz(1)-1,min_ind_kz2:min_ind_kz2+sz(2)-1);
image_mask = ifft2c(k_mask);
% Threshold to binarise re-scaled image_mask
image_mask(abs(image_mask) < prctile(abs(image_mask),(1-sum(image_mask_full,'all')./(size(image_mask_full,1)*size(image_mask_full,2)))*1e2,'all')) = 0;
image_mask(abs(image_mask) >= prctile(abs(image_mask),(1-sum(image_mask_full,'all')./(size(image_mask_full,1)*size(image_mask_full,2)))*1e2,'all')) = 1;
%figure(); imagesc(abs(image_mask)); axis image off; title('Mask')

%% Perform image recon for accelerated and unaccelerated images
tic
% Calculate extent of ESPIRiT calibration region
centre = [sx,sy]/2;
calib1 = centre(1)-floor((calibsize(1)-1)/2):centre(1)+ceil((calibsize(1)-1)/2);
calib2 = centre(2)-floor((calibsize(2)-1)/2):centre(2)+ceil((calibsize(2)-1)/2);

% Calculate sensitivty maps
rx_sens = tx_espirit(permute(k_Relativec(calib1,calib2,:,:),[1,2,4,3]), sz, kernel, eig_thresh);
tx_sens = tx_espirit(k_Relativec(calib1,calib2,:,:), sz, kernel, eig_thresh);


if (sz(1) ~= sx || sz(2) ~= sy) && any(size(k_Relativec,[1,2]) ~= sz)
    %Hanning filter and zero-pad k-space
    Hanning2D = hann(sx)*hann(sy)';
    k_Relativec = padarray(padarray(k_Relativec.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
    k_Ref_Shims_mTxc = padarray(padarray(k_Ref_Shims_mTxc.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
    k_Prep_Shims_mTxc = padarray(padarray(k_Prep_Shims_mTxc.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
end

%[rx_sens2,~,~] = ESPIRiTmaps(sum(k_Relativec(calib1,calib2,:,:),4),kernel,eig_thresh,eig_thresh2,0); % Sum across transmit dimension
%[tx_sens,~,~] = ESPIRiTmaps(sum(k_Relativec,3),kernel,eig_thresh,eig_thresh2,0); % Sum across receive dimension
%figure(); imagesc(imtile(abs(squeeze(rx_sens))),[0 1])
%figure(); imagesc(imtile(abs(squeeze(rx_sens2(:,:,:,2)))),[0 1])

% FFT back to image space
im_Relativec = ifft2c(k_Relativec);
im_Ref_Shims_mTxc = ifft2c(k_Ref_Shims_mTxc);
im_Prep_Shims_mTxc = ifft2c(k_Prep_Shims_mTxc);

% Combine sensitivity maps with images, sum across receive
% channels without acceleration
Relative_Images = squeeze(sum(bsxfun(@times,im_Relativec,conj(rx_sens)),3));
Ref_Shims = squeeze(sum(bsxfun(@times,im_Ref_Shims_mTxc,conj(rx_sens)),3));
Prep_Shims = squeeze(sum(bsxfun(@times,im_Prep_Shims_mTxc,conj(rx_sens)),3));

plotsupportingfigureS1(Relative_Images,Ref_Shims,Prep_Shims) % Supporting Figure 1

for iter_n  = 1:size(niters_array,2)
    niters = niters_array(iter_n);
    load([Folder1,filename,'.mat']); % load in kspace data
    filename2 = [filename,'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat']; % filename for reconstructed image data
    if exist([Folder2,filename2],'file') == 2
        disp('Data already exists. Loading in')
        load([Folder2,filename2]);
    else
        % pre-allocate
        Rel_Maps = zeros([sz, 8]);
        Rel_Maps_acc = zeros([sz, 8]);
        Maps = zeros([sz, 8, NRepeats,3,size(accelerations,2)]);
        Maps_acc = zeros([sz, 8, NRepeats,3,size(accelerations,2)]);
        rx_sens_acc = zeros([sz, 8, NRepeats,3,size(accelerations,2)]);
        tx_sens_acc = zeros([sz, 8, NRepeats,3,size(accelerations,2)]);
        
        % Calculate sensitivity maps from accelerated relative k-space
        parfor accel_ind = 1:size(accelerations,2)
            for same_masks = 2%1:3 % Code used for testing different alternative masking strategies
                for repeat_n = 1:NRepeats
                    % Calculate sensitivity maps from cropped k-space for
                    % accelerated data only from Relative images
                    rx_sens_acc(:,:,:,repeat_n,same_masks,accel_ind) = tx_espirit(permute(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind),[1,2,4,3]), sz, kernel, eig_thresh);
                    tx_sens_acc(:,:,:,repeat_n,same_masks,accel_ind) = tx_espirit(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind), sz, kernel, eig_thresh);
                    %[rx_sens_acc(:,:,:,:,repeat_n,same_masks,accel_ind),~,cgESPIRiTrecon(:,:,:,:,repeat_n,same_masks,accel_ind)] = ESPIRiTmaps(sum(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind),4),kernel,eig_thresh,eig_thresh2,0,cat(4,refUS(:,:,:,:,repeat_n,same_masks,accel_ind),prepUS(:,:,:,:,repeat_n,same_masks,accel_ind))); % Sum across transmit dimension
                    %[tx_sens_acc(:,:,:,:,repeat_n,same_masks,accel_ind),~,~] = ESPIRiTmaps(sum(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind),3),kernel,eig_thresh,eig_thresh2,0); % Sum across receive dimension
                end
            end
        end
        
        if (sz(1) ~= sx || sz(2) ~= sy) && any(size(recon_rel,[1,2]) ~= sz)
            % Hanning filter to padd k-space
            Hanning2D = hann(sx)*hann(sy)';
            recon_rel = padarray(padarray(recon_rel.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
            if strcmp(Recon_Type,'ESPIRiT') ||  strcmp(Recon_Type,'GRAPPA')
                refUS = padarray(padarray(refUS.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
                prepUS = padarray(padarray(prepUS.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
                recon_ref = zeros([sz,8,size(k_Ref_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]); % pre-allocate
                recon_prep = zeros([sz,8,size(k_Prep_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]); % pre-allocate
            else
                recon_ref = padarray(padarray(recon_ref.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
                recon_prep = padarray(padarray(recon_prep.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
            end
        end
        
        % Perform cgESPIRiT reconstruction on reference and prepared
        for accel_ind = 1:size(accelerations,2)
            for same_masks = 2%1:3 % Code used for testing different alternative masking strategies
                for repeat_n = 1:NRepeats
                    if strcmp(Recon_Type,'ESPIRiT')
                        % ESPIRiT CG reconstruction of reference/prepared images
                        cgiters = 12;
                        SNS = ESPIRiT(rx_sens_acc(:,:,:,repeat_n,same_masks,accel_ind));
                        [recon_ref(:,:,1,repeat_n,same_masks,accel_ind),~]  = cgESPIRiT(refUS(:,:,:,1,repeat_n,same_masks,accel_ind), SNS, cgiters,0.01 ,refUS(:,:,:,1,repeat_n,same_masks,accel_ind)*0);
                        [recon_ref(:,:,2,repeat_n,same_masks,accel_ind),~]  = cgESPIRiT(refUS(:,:,:,2,repeat_n,same_masks,accel_ind), SNS, cgiters,0.01 ,refUS(:,:,:,2,repeat_n,same_masks,accel_ind)*0);
                        [recon_prep(:,:,1,repeat_n,same_masks,accel_ind),~]  = cgESPIRiT(prepUS(:,:,:,1,repeat_n,same_masks,accel_ind), SNS, cgiters,0.01 ,prepUS(:,:,:,1,repeat_n,same_masks,accel_ind)*0);
                        [recon_prep(:,:,2,repeat_n,same_masks,accel_ind),~]  = cgESPIRiT(prepUS(:,:,:,2,repeat_n,same_masks,accel_ind), SNS, cgiters,0.01 ,prepUS(:,:,:,2,repeat_n,same_masks,accel_ind)*0);
                    elseif strcmp(Recon_Type,'GRAPPA')
                        % GRAPPA reconstruction of reference/prepared images
                        calibregion = sum(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind),4);
                        recon_ref(:,:,:,1,repeat_n,same_masks,accel_ind)  =  GRAPPA(refUS(:,:,:,1,repeat_n,same_masks,accel_ind),calibregion,[5,5],0.01);
                        recon_ref(:,:,:,2,repeat_n,same_masks,accel_ind)  =  GRAPPA(refUS(:,:,:,2,repeat_n,same_masks,accel_ind),calibregion,[5,5],0.01);
                        recon_prep(:,:,:,1,repeat_n,same_masks,accel_ind)  =  GRAPPA(prepUS(:,:,:,1,repeat_n,same_masks,accel_ind),calibregion,[5,5],0.01);
                        recon_prep(:,:,:,2,repeat_n,same_masks,accel_ind)  =  GRAPPA(prepUS(:,:,:,2,repeat_n,same_masks,accel_ind),calibregion,[5,5],0.01);
                    end
                end
            end
        end
        
        % FFT back to image space
        im_Rel_acc = ifft2c(recon_rel);
        im_Ref_acc = ifft2c(recon_ref);
        im_Prep_acc = ifft2c(recon_prep);
        
        % Combine sensitivity maps with images, sum across receive
        % channels with acceleration
        Relative_Images_acc = squeeze(sum(bsxfun(@times,im_Rel_acc,conj(permute(rx_sens_acc,[1 2 3 7 4 5 6]))),3));
        Ref_Shims_acc = squeeze(sum(bsxfun(@times,im_Ref_acc,conj(permute(rx_sens_acc,[1 2 3 7 4 5 6]))),3));
        Prep_Shims_acc = squeeze(sum(bsxfun(@times,im_Prep_acc,conj(permute(rx_sens_acc,[1 2 3 7 4 5 6]))),3));
        
        %figure('color','w'); imagesc(imtile(abs(squeeze(Relative_Images_acc(:,:,8,1,2,4))))); axis image off
        %figure('color','w'); imagesc(imtile(abs(squeeze(Ref_Shims_acc(:,:,1,1,2,3))))); axis image off
        %figure('color','w'); imagesc(imtile(abs(squeeze(Prep_Shims_acc(:,:,1,1,2,3))))); axis image off
        
        %nrmse(Relative_Images_acc(:,:,:,1,2,4),Relative_Images_acc(:,:,:,1,2,1))
        %nrmse(Ref_Shims_acc(:,:,:,1,2,4),Ref_Shims_acc(:,:,:,1,2,1))
        %nrmse(Prep_Shims_acc(:,:,:,1,2,4),Prep_Shims_acc(:,:,:,1,2,1))
        
        % Calculate relative maps from transmit sensitivities
        % or relative images
        Rel_Maps = squeeze(tx_sens./sum(abs(tx_sens),3));
        Rel_Maps_acc = squeeze(tx_sens_acc./sum(abs(tx_sens_acc),3));
        %Rel_Maps = squeeze(Relative_Images./sum(abs(Relative_Images),3));
        %Rel_Maps_acc = squeeze(Relative_Images_acc./sum(abs(Relative_Images_acc),3));
        
        % Absolute shim maps (un-accelerated)
        FA_Shim1 = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims(:,:,1:2:end) ./ Ref_Shims(:,:,1:2:end)))));
        FA_Shim2 = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims(:,:,2:2:end) ./ Ref_Shims(:,:,2:2:end)))));
        % Absolute shim maps (accelerated)
        FA_Shim1_acc = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims_acc(:,:,1,:,:,:) ./ Ref_Shims_acc(:,:,1,:,:,:)))));
        FA_Shim2_acc = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims_acc(:,:,2,:,:,:) ./ Ref_Shims_acc(:,:,2,:,:,:)))));
        %figure('color','w'); imagesc(imtile({abs(FA_Shim1),abs(FA_Shim2)}),[0 150]); axis image off
        %figure('color','w'); imagesc(imtile({abs(FA_Shim1_acc(:,:,1,2,1)),abs(FA_Shim2_acc(:,:,1,2,1))}),[0 150]); axis image off
        %figure('color','w'); imagesc(imtile({abs(FA_Shim1),abs(FA_Shim2),abs(FA_Shim1_acc(:,:,1,2,3)),abs(FA_Shim2_acc(:,:,1,2,3))}),[0 150]); axis image off; title('Fully Sampled vs Under Sampled');
        
        % Now calculate absolute single channel maps (un-accelerated)
        Rel_Shim1 = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps,permute(Shim_Setting1,[1 3 2])),3));
        Rel_Shim2 = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps,permute(Shim_Setting2,[1 3 2])),3));
        % Now calculate absolute single channel maps (accelerated)
        Rel_Shim1_acc = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps_acc,permute(Shim_Setting1,[1 3 2])),3));
        Rel_Shim2_acc = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps_acc,permute(Shim_Setting2,[1 3 2])),3));
        %figure(); imagesc(imtile({abs(Rel_Shim1),abs(Rel_Shim2)}))
        
        % Divide FA shims by relative shims to get absolute scaling factors
        % (un-accelerated)
        c1 = abs(FA_Shim1./Rel_Shim1); % figure(); imagesc(imtile({abs(FA_Shim1),abs(Rel_Shim1)})) % Units degrees
        c2 = abs(FA_Shim2./Rel_Shim2); % figure(); imagesc(imtile({abs(FA_Shim2),abs(Rel_Shim2)})) % Units degrees
        % Divide FA shims by relative shims to get absolute scaling factors
        % (accelerated)
        c1_acc = abs(FA_Shim1_acc./Rel_Shim1_acc);
        c2_acc = abs(FA_Shim2_acc./Rel_Shim2_acc);
        %figure(); imagesc(imtile({abs(c1),abs(c2)}))
        
        % Weighted combination of scaling factors
        m = 4;
        C = ((c1.*abs(Rel_Shim1).^m) + (c2.*abs(Rel_Shim2).^m))./ (abs(Rel_Shim1).^m + abs(Rel_Shim2).^m);
        C_acc = ((c1_acc.*abs(Rel_Shim1_acc).^m) + (c2_acc.*abs(Rel_Shim2_acc).^m))./ (abs(Rel_Shim1_acc).^m + abs(Rel_Shim2_acc).^m); % Note B1TIAMO equation is incorrect in paper? m inside brackets
        C(isnan(C)) = 0; C_acc(isnan(C_acc)) = 0; % Remove NaNs
        
        % Scale relative maps to absolute maps
        Maps = image_mask.*C.*Rel_Maps;
        Maps_acc = image_mask.*repmat(permute(C_acc,[1 2 6 3 4 5]),[1 1 8]).*Rel_Maps_acc;
        
        % Remove nans
        Maps(isnan(Maps))=0;
        Maps_acc(isnan(Maps_acc))=0;
        
        % There may be different masking for the fully sampled vs accelerated.
        %Hence, here we find a mask where both have values and apply this.
        mask1 = zeros(size(Maps)); mask2 = mask1; mask3 = mask1;
        mask1(Maps~=0) =1; mask2(Maps_acc(:,:,:,1,2,1)~=0) =1;
        mask3(mask1&mask2) = 1;
        Maps = mask3.*Maps; Maps_acc = mask3.*Maps_acc;
        
        save([Folder2,filename2],'Maps','Maps_acc');
    end
    disp(['Iteration ',num2str(niters_array(iter_n)), ' completed for all acceleration factors ', datestr(now, 'dd/mm/yy-HH:MM')])
end



toc
disp('Finished Simulations.')
end