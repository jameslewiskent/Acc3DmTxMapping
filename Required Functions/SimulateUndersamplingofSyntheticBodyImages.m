function [] = SimulateUndersamplingofSyntheticBodyImages(settings)
load('SyntheticBodyImages.mat','Relative_Rx','Ref_Shims','Prep_Shims','image_mask_full','image_heartmask_full');

sx = settings.sx;
sy = settings.sy;
sz = settings.sz;
calib = settings.calib;
niters_array = settings.niters_array;
accelerations = settings.accelerations;
[Enc_Mat,~] = CalcEncMat(settings.Enc_Scheme);
Shim_Setting1 = Enc_Mat(9,:); % CP
Shim_Setting2 = Enc_Mat(10,:); % CP2
NRepeats = settings.NRepeats;
pSNR = settings.pSNR;

[masks] = genCSmasks(settings); % generate undersampling masks if neccesary

kernel = [5 5]; % kernel size
eig_thresh = 0.02; % Threshold for picking singular vectors of the calibration matrix (relative to largest singular value)

location = ['Data',filesep,'Synthetic Body Simulation Results',filesep];
Folder1 = [location,'AcceleratedData',filesep];
Folder2 = [location,'ReconData',filesep];
if exist(Folder1,'dir') ~= 7 || exist(Folder2,'dir') ~= 7
    mkdir(Folder1); mkdir(Folder2);
end

%% Crop out arms (optional)
%figure(); imagesc(imtile(image_mask_full))
image_mask_full(:,1:8,:,:)= 0; image_mask_full(:,170:178,:,:)= 0;
%figure(); imagesc(imtile(image_mask_full))
Relative_Rx = Relative_Rx.*image_mask_full;
Ref_Shims = Ref_Shims.*image_mask_full;
Prep_Shims = Prep_Shims.*image_mask_full;

%% FFT and crop k-space to desired size (sx, sy)
kdata_Relative = fft2c(Relative_Rx);
kdata_Ref_Shims = fft2c(Ref_Shims);
kdata_Prep_Shims = fft2c(Prep_Shims);

% Crop kspace to sz
min_ind_kx = ceil(round((size(kdata_Relative,1)+1)/2) - sz(1)/2);
min_ind_ky = ceil(round((size(kdata_Relative,2)+1)/2) - sz(2)/2);
k_Relativec = kdata_Relative(min_ind_kx:min_ind_kx+sz(1)-1,min_ind_ky:min_ind_ky+sz(2)-1,:,:);
k_Ref_Shimsc = kdata_Ref_Shims(min_ind_kx:min_ind_kx+sz(1)-1,min_ind_ky:min_ind_ky+sz(2)-1,:,:);
k_Prep_Shimsc = kdata_Prep_Shims(min_ind_kx:min_ind_kx+sz(1)-1,min_ind_ky:min_ind_ky+sz(2)-1,:,:);
clear min_ind_kx min_ind_ky

% Use to crop kspace to [sx,sy]
min_ind_kx = ceil(round((size(k_Relativec,1)+1)/2) - sx/2);
min_ind_ky = ceil(round((size(k_Relativec,2)+1)/2) - sy/2);
crop_region1 = min_ind_kx:min_ind_kx+sx-1;
crop_region2 = min_ind_ky:min_ind_ky+sy-1;

% define noise levels
std1 = max(abs(k_Relativec),[],'all')./db2mag(pSNR); std2 = max(abs(k_Ref_Shimsc),[],'all')./db2mag(pSNR); std3 = max(abs(k_Prep_Shimsc),[],'all')./db2mag(pSNR);
%% Under-sample images using under-sampling masks, then reconstruct k-space using admm
tic
for iter_n  = 1:size(niters_array,2)
    niters = niters_array(iter_n);
    filename = ['Simulated_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats)];
    
    % pre-allocate
    recon_rel = zeros([sx,sy,8,size(k_Relativec,4),NRepeats,3,size(accelerations,2)]);
    recon_ref = zeros([sx,sy,8,size(k_Ref_Shimsc,4),NRepeats,3,size(accelerations,2)]);
    recon_prep = zeros([sx,sy,8,size(k_Prep_Shimsc,4),NRepeats,3,size(accelerations,2)]);
    
    if exist([Folder1,filename,'.mat'],'file') == 2
        disp(['File ', filename,'.mat', ' already exists, loading in data'])
        load([Folder1,filename,'.mat'],'recon_rel','recon_ref','recon_prep');
    else
        for accel_ind = 1:size(accelerations,2)
            kdata_Relativec_acc = zeros(sx,sy,8,8);
            kdata_Ref_Shimsc_acc = zeros(sx,sy,8,2);
            kdata_Prep_Shimsc_acc = zeros(sx,sy,8,2);
            for same_masks = 2 %1:3 Code used for testing different alternative masking strategies
                for repeat_n = 1:NRepeats
                    % Corrupt with noise
                    k_Relativecn = k_Relativec(crop_region1,crop_region2,:,:) + std1*randn(size(k_Relativec(crop_region1,crop_region2,:,:))) + 1i*std1*randn(size(k_Relativec(crop_region1,crop_region2,:,:)));
                    k_Ref_Shimscn = k_Ref_Shimsc(crop_region1,crop_region2,:,:) + std2*randn(size(k_Ref_Shimsc(crop_region1,crop_region2,:,:))) + 1i*std2*randn(size(k_Ref_Shimsc(crop_region1,crop_region2,:,:)));
                    k_Prep_Shimscn  = k_Prep_Shimsc(crop_region1,crop_region2,:,:)  + std3*randn(size(k_Prep_Shimsc(crop_region1,crop_region2,:,:))) + 1i*std3*randn(size(k_Prep_Shimsc(crop_region1,crop_region2,:,:)));
                    
                    % Choose masks
                    if same_masks == 1 % all Tx modes, all ref/prep images use same masks
                        rel_mask_n = repeat_n;
                        ref_mask_n = repeat_n;
                        prep_mask_n = repeat_n;
                    elseif same_masks == 2 % Tx modes use different masks, ref/prep images use same masks (paired)
                        n_refprep_masks = size(k_Prep_Shimsc,4); n_rel_masks = size(k_Relativec,4); total_per_rep = n_refprep_masks + n_rel_masks;
                        rep_base = ((repeat_n-1)*total_per_rep)+1;
                        rel_mask_n = (rep_base:rep_base+n_rel_masks-1);
                        ref_mask_n = (rep_base+n_rel_masks:rep_base+n_rel_masks+n_refprep_masks-1);
                        prep_mask_n = ref_mask_n; % reference image uses same mask as prepared
                    elseif same_masks == 3 % all Tx modes, all ref/prep images use different masks
                        n_refprep_masks = size(k_Prep_Shimsc,4); n_rel_masks = size(k_Relativec,4); total_per_rep = 2*n_refprep_masks + n_rel_masks;
                        rep_base = ((repeat_n-1)*total_per_rep)+1;
                        rel_mask_n = (rep_base:rep_base+n_rel_masks-1);
                        ref_mask_n = (rep_base+n_rel_masks:rep_base+n_rel_masks+n_refprep_masks-1);
                        prep_mask_n = (rep_base+n_rel_masks+n_refprep_masks:rep_base+n_rel_masks+2*n_refprep_masks-1);
                    end
                    % Apply masks
                    kdata_Relativec_acc = bsxfun(@times,k_Relativecn,masks(:,:,accel_ind,rel_mask_n));
                    kdata_Ref_Shimsc_acc = bsxfun(@times,k_Ref_Shimscn,masks(:,:,accel_ind,ref_mask_n));
                    kdata_Prep_Shimsc_acc = bsxfun(@times,k_Prep_Shimscn,masks(:,:,accel_ind,prep_mask_n));
                    
                    % TxLR recon
                    joint_recon = admm_txlr(double(cat(4,kdata_Relativec_acc,kdata_Ref_Shimsc_acc,kdata_Prep_Shimsc_acc)), kernel, niters, [50 50]);
                    recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind) = joint_recon(:,:,:,1:8);
                    recon_ref(:,:,:,:,repeat_n,same_masks,accel_ind) = joint_recon(:,:,:,9:10);
                    recon_prep(:,:,:,:,repeat_n,same_masks,accel_ind) = joint_recon(:,:,:,11:12);
                    
                end
            end
            disp(['Acceleration Factor of ',num2str(accelerations(accel_ind)), ' completed ', datestr(now, 'dd/mm/yy-HH:MM')])
        end
        save([Folder1,filename,'.mat'],'recon_rel','recon_ref','recon_prep')
    end
end
toc

clear k_Relativecn k_Ref_Shimscn k_Prep_Shimscn
%% Define mask
image_mask = ScaleMask(image_mask_full,sz);
heart_mask = ScaleMask(image_heartmask_full,sz);

%% Perform image recon for accelerated and unaccelerated images
tic

TukeyRot = RotArray(tukeywin(sx+2,0.7)); TukeyRot = TukeyRot(2:end-1,2:end-1); % Tukey Filter

% Calculate sensitivty maps
rx_sens = tx_espirit(permute(k_Relativec,[1,2,4,3]), sz, kernel, eig_thresh);
tx_sens = tx_espirit(k_Relativec, sz, kernel, eig_thresh);

if (sz(1) ~= sx || sz(2) ~= sy) && any(size(k_Relativec,[1,2]) ~= sz)
    k_Relativec = padarray(padarray(k_Relativec.*TukeyRot,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
    k_Ref_Shimsc = padarray(padarray(k_Ref_Shimsc.*TukeyRot,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
    k_Prep_Shimsc = padarray(padarray(k_Prep_Shimsc.*TukeyRot,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
end

% FFT back to image space
im_Relativec = ifft2c(k_Relativec);
im_Ref_Shimsc = ifft2c(k_Ref_Shimsc);
im_Prep_Shimsc = ifft2c(k_Prep_Shimsc);

% Combine sensitivity maps with images, sum across receive
% channels without acceleration
Relative_Images = squeeze(sum(bsxfun(@times,im_Relativec,conj(rx_sens)),3));
Ref_Shims = squeeze(sum(bsxfun(@times,im_Ref_Shimsc,conj(rx_sens)),3));
Prep_Shims = squeeze(sum(bsxfun(@times,im_Prep_Shimsc,conj(rx_sens)),3));

plotsupportingfigureS1(Relative_Images,Ref_Shims,Prep_Shims) % Supporting Figure 1

for iter_n  = 1:size(niters_array,2)
    niters = niters_array(iter_n);
    filename = ['Simulated_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats)];
    load([Folder1,filename,'.mat'],'recon_rel','recon_ref','recon_prep'); % load in kspace data
    filename2 = [filename,'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat']; % filename for reconstructed image data
    if exist([Folder2,filename2],'file') == 2
        disp('Data already exists. Loading in')
        load([Folder2,filename2],'Maps','Maps_acc');
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
                    rx_sens_acc(:,:,:,repeat_n,same_masks,accel_ind) = tx_espirit(permute(recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind),[1,2,4,3]), sz, kernel, eig_thresh);
                    tx_sens_acc(:,:,:,repeat_n,same_masks,accel_ind) = tx_espirit(recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind), sz, kernel, eig_thresh);
                end
            end
        end
        
        if (sz(1) ~= sx || sz(2) ~= sy) && any(size(recon_rel,[1,2]) ~= sz)
            recon_rel = padarray(padarray(recon_rel.*TukeyRot,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
            recon_ref = padarray(padarray(recon_ref.*TukeyRot,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
            recon_prep = padarray(padarray(recon_prep.*TukeyRot,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
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
                
        % Now calculate absolute single channel maps (un-accelerated)
        Rel_Shim1 = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps,permute(Shim_Setting1,[1 3 2])),3));
        Rel_Shim2 = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps,permute(Shim_Setting2,[1 3 2])),3));
        % Now calculate absolute single channel maps (accelerated)
        Rel_Shim1_acc = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps_acc,permute(Shim_Setting1,[1 3 2])),3));
        Rel_Shim2_acc = squeeze(image_mask.*sum(bsxfun(@times,Rel_Maps_acc,permute(Shim_Setting2,[1 3 2])),3));
        
        % Divide FA shims by relative shims to get absolute scaling factors
        % (un-accelerated)
        c1 = abs(FA_Shim1./Rel_Shim1); % figure(); imagesc(imtile({abs(FA_Shim1),abs(Rel_Shim1)})) % Units degrees
        c2 = abs(FA_Shim2./Rel_Shim2); % figure(); imagesc(imtile({abs(FA_Shim2),abs(Rel_Shim2)})) % Units degrees
        % Divide FA shims by relative shims to get absolute scaling factors
        % (accelerated)
        c1_acc = abs(FA_Shim1_acc./Rel_Shim1_acc);
        c2_acc = abs(FA_Shim2_acc./Rel_Shim2_acc);
        
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
        % Hence, here we find a mask where both have values and apply this.
        mask1 = zeros(size(Maps)); mask2 = mask1; mask3 = mask1;
        mask1(Maps~=0) =1; mask2(Maps_acc(:,:,:,1,2,1)~=0) =1;
        mask3(mask1&mask2) = 1;
        Maps = mask3.*Maps; Maps_acc = mask3.*Maps_acc;
        
        save([Folder2,filename2],'Maps','Maps_acc','heart_mask');
    end
    disp(['Iteration ',num2str(niters_array(iter_n)), ' completed for all acceleration factors ', datestr(now, 'dd/mm/yy-HH:MM')])
end

toc
disp('Finished Simulations.')
end