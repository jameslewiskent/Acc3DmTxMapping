function [] = SimulateUndersamplingofSyntheticBodyImages(sx,sy,sz,calib,niters,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

load('SyntheticBodyImages.mat','Relative_Rx','Ref_Shims_mTx','Prep_Shims_mTx','image_mask_full');

kernel = [5,5]; % ESPIRiT kernel
eig_thresh = 0.02; % ESPIRiT eigenvalue threshold
calibsize = [sx,sy]; % ESPIRiT calibration size

location = ['Data',filesep,'Synthetic Body Simulation Results',filesep];
Folder1 = [location,'AcceleratedData',filesep];
Folder2 = [location,'ReconData',filesep];
%% Crop out arms
% %figure(); imagesc(imtile(image_mask_full))
% image_mask_full(:,1:8,:,:)= 0; image_mask_full(:,170:178,:,:)= 0; 
% %figure(); imagesc(imtile(image_mask_full))
% Relative_Rx = Relative_Rx.*image_mask_full;
% Ref_Shims_mTx = Ref_Shims_mTx.*image_mask_full;
% Prep_Shims_mTx = Prep_Shims_mTx.*image_mask_full;

%% FFT and crop k-space to desired size (sx, sy)
kdata_Relative = fftshift(fftshift(fft(fft(ifftshift(ifftshift(Relative_Rx,1),2),[],1),[],2),1),2);
kdata_Ref_Shims_mTx = fftshift(fftshift(fft(fft(ifftshift(ifftshift(Ref_Shims_mTx,1),2),[],1),[],2),1),2);
kdata_Prep_Shims_mTx = fftshift(fftshift(fft(fft(ifftshift(ifftshift(Prep_Shims_mTx,1),2),[],1),[],2),1),2);

% Corrupt with noise
std1 = max(abs(kdata_Relative),[],'all')./db2mag(pSNR); std2 = max(abs(kdata_Ref_Shims_mTx),[],'all')./db2mag(pSNR); std3 = max(abs(kdata_Prep_Shims_mTx),[],'all')./db2mag(pSNR);
kdata_Relative = kdata_Relative + std1*randn(size(kdata_Relative)) + 1i*std1*randn(size(kdata_Relative));
kdata_Ref_Shims_mTx = kdata_Ref_Shims_mTx + std2*randn(size(kdata_Ref_Shims_mTx)) + 1i*std2*randn(size(kdata_Ref_Shims_mTx));
kdata_Prep_Shims_mTx  = kdata_Prep_Shims_mTx  + std3*randn(size(kdata_Prep_Shims_mTx)) + 1i*std3*randn(size(kdata_Prep_Shims_mTx));

% Crop kspace to sx * sy
min_ind_kx = ceil(round((size(kdata_Relative,1)+1)/2) - sx/2);
min_ind_ky = ceil(round((size(kdata_Relative,2)+1)/2) - sy/2);
k_Relativec = kdata_Relative(min_ind_kx:min_ind_kx+sx-1,min_ind_ky:min_ind_ky+sy-1,:,:);
k_Ref_Shims_mTxc = kdata_Ref_Shims_mTx(min_ind_kx:min_ind_kx+sx-1,min_ind_ky:min_ind_ky+sy-1,:,:);
k_Prep_Shims_mTxc = kdata_Prep_Shims_mTx(min_ind_kx:min_ind_kx+sx-1,min_ind_ky:min_ind_ky+sy-1,:,:);

%% Under-sample images using under-sampling masks, then reconstruct k-space using admm
tic
for iter_n  = 1:size(niters,2)
    niters = niters(iter_n);
    filename = ['Simulated_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats)];
    
    % pre-allocate
    recon_rel = zeros([sx,sy,8,size(k_Relativec,4),NRepeats,3,size(accelerations,2)]);
    recon_ref = zeros([sx,sy,8,size(k_Ref_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]);
    recon_prep = zeros([sx,sy,8,size(k_Prep_Shims_mTxc,4),NRepeats,3,size(accelerations,2)]);
    err_kdata_rel = zeros([NRepeats,3,size(accelerations,2)]);
    err_kdata_ref = zeros([NRepeats,3,size(accelerations,2)]);
    err_kdata_prep = zeros([NRepeats,3,size(accelerations,2)]);
    
    if exist([Folder1,filename,'.mat'],'file') == 2
        disp(['File ', filename,'.mat', ' already exists, loading in data'])
        load([Folder1,filename,'.mat']);
    else
        parfor accel_ind = 1:size(accelerations,2)
            for same_masks = 2 %1:3 Code used for testing different alternative masking strategies
                for repeat_n = 1:NRepeats
                    
                    % Controls same/different masks for Tx modes and
                    % ref/prepared images
                    if same_masks == 1 % all Tx modes, all ref/prep images use same masks
                        rel_mask_n = repeat_n;
                        ref_mask_n = repeat_n;
                        prep_mask_n = repeat_n;
                        
                        kdata_Relativec_acc = bsxfun(@times,k_Relativec,masks(:,:,accel_ind,rel_mask_n));
                        kdata_Ref_Shims_mTxc_acc = bsxfun(@times,k_Ref_Shims_mTxc,masks(:,:,accel_ind,ref_mask_n));
                        kdata_Prep_Shims_mTxc_acc = bsxfun(@times,k_Prep_Shims_mTxc,masks(:,:,accel_ind,prep_mask_n));
                    elseif same_masks == 2 % Tx modes use different masks, ref/prep images use same masks (paired)
                        n_refprep_masks = size(k_Prep_Shims_mTxc,4); n_rel_masks = size(k_Relativec,4); count = 1; total_per_rep = n_refprep_masks + n_rel_masks;
                        rep_base = ((repeat_n-1)*total_per_rep)+1;
                        rel_mask_n = (rep_base:rep_base+n_rel_masks-1);
                        ref_mask_n = (rep_base+n_rel_masks:rep_base+n_rel_masks+n_refprep_masks-1);
                        prep_mask_n = ref_mask_n; % reference image uses same mask as prepared
                        
                        kdata_Relativec_acc = bsxfun(@times,k_Relativec,masks(:,:,accel_ind,rel_mask_n));
                        kdata_Ref_Shims_mTxc_acc = bsxfun(@times,k_Ref_Shims_mTxc,masks(:,:,accel_ind,ref_mask_n));
                        kdata_Prep_Shims_mTxc_acc = bsxfun(@times,k_Prep_Shims_mTxc,masks(:,:,accel_ind,prep_mask_n));
                    elseif same_masks == 3 % all Tx modes, all ref/prep images use different masks
                        n_refprep_masks = size(k_Prep_Shims_mTxc,4); n_rel_masks = size(k_Relativec,4); count = 1; total_per_rep = 2*n_refprep_masks + n_rel_masks;
                        rep_base = ((repeat_n-1)*total_per_rep)+1;
                        rel_mask_n = (rep_base:rep_base+n_rel_masks-1);
                        ref_mask_n = (rep_base+n_rel_masks:rep_base+n_rel_masks+n_refprep_masks-1);
                        prep_mask_n = (rep_base+n_rel_masks+n_refprep_masks:rep_base+n_rel_masks+2*n_refprep_masks-1);
                        
                        kdata_Relativec_acc = bsxfun(@times,k_Relativec,masks(:,:,accel_ind,rel_mask_n));
                        kdata_Ref_Shims_mTxc_acc = bsxfun(@times,k_Ref_Shims_mTxc,masks(:,:,accel_ind,ref_mask_n));
                        kdata_Prep_Shims_mTxc_acc = bsxfun(@times,k_Prep_Shims_mTxc,masks(:,:,accel_ind,prep_mask_n));
                    end
                    
                    
                    % Joint recon of rel, reference & prep images by concatenating along Tx dimension
                    %joint_data = cat(4,kdata_Relativec_acc,kdata_Ref_Shims_mTxc_acc,kdata_Prep_Shims_mTxc_acc);
                    %recon_joint = admm_txlr(double(joint_data), kernel, niters, [50 50]);
                    
                    recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind) = admm_txlr(double(kdata_Relativec_acc), kernel, niters, [50 50]);

                    %recon_rel(:,:,:,:,repeat_n,same_masks,accel_ind) = recon_joint(:,:,:,1:8);
                    %recon_ref(:,:,:,:,repeat_n,same_masks,accel_ind) = recon_joint(:,:,:,9:9+size(k_Ref_Shims_mTxc,4)-1);
                    %recon_prep(:,:,:,:,repeat_n,same_masks,accel_ind) = recon_joint(:,:,:,9+size(k_Prep_Shims_mTxc,4):end);
               end
            end
            disp(['Acceleration Factor of ',num2str(accelerations(accel_ind)), ' completed ', datestr(now, 'dd/mm/yy-HH:MM')])
        end
        save([Folder1,filename,'.mat'],'err_kdata_rel','recon_rel','recon_ref','recon_prep','err_kdata_ref','err_kdata_prep')
    end
end
toc
%% Define mask
% reduce size of image mask
image_mask = fftshift(fftshift(fft2(ifftshift(ifftshift(image_mask_full,1),2)),1),2);
min_ind_kz1 = round(round((size(image_mask,1)+1)/2) - sz(1)/2);
min_ind_kz2 = round(round((size(image_mask,2)+1)/2) - sz(2)/2);
image_mask = image_mask(min_ind_kz1:min_ind_kz1+sz(1)-1,min_ind_kz2:min_ind_kz2+sz(2)-1);
image_mask = fftshift(fftshift(ifft2(ifftshift(ifftshift(image_mask,1),2)),1),2);
% Threshold to binarise re-scaled image_mask
image_mask(abs(image_mask) < prctile(abs(image_mask),(1-sum(image_mask_full,'all')./(size(image_mask_full,1)*size(image_mask_full,2)))*1e2,'all')) = 0;
image_mask(abs(image_mask) >= prctile(abs(image_mask),(1-sum(image_mask_full,'all')./(size(image_mask_full,1)*size(image_mask_full,2)))*1e2,'all')) = 1;
%figure(); imagesc(abs(image_mask)); axis image off; title('Mask')

%% Perform image recon for accelerated and unaccelerated images
tic
% Calculate extent of ESPIRiT calibration region
centre = [sx, sy]/2;
calib1 = centre(1)-floor((calibsize(1)-1)/2):centre(1)+ceil((calibsize(1)-1)/2);
calib2 = centre(2)-floor((calibsize(2)-1)/2):centre(2)+ceil((calibsize(2)-1)/2);

% Calculate sensitivity maps from cropped relative map k-space (for
% specific calibration region)
rx_sens = tx_espirit(permute(k_Relativec(calib1,calib2,:,:),[1,2,4,3]), sz, kernel, eig_thresh);
tx_sens = tx_espirit(k_Relativec(calib1,calib2,:,:), sz, kernel, eig_thresh);

%figure(); imagesc(imtile(abs(squeeze(rx_sens(:,:,:,1)))),[0 1])
%figure(); imagesc(imtile(abs(squeeze(rx_sens(:,:,:,2)))),[0 1])

if (sz(1) ~= sx || sz(2) ~= sy) && any(size(k_Relativec,[1,2]) ~= sz)
    % Hanning filter and zero-pad k-space
    Hanning2D = hann(sx)*hann(sy)';
    k_Relativec = padarray(padarray(k_Relativec.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
    k_Ref_Shims_mTxc = padarray(padarray(k_Ref_Shims_mTxc.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
    k_Prep_Shims_mTxc = padarray(padarray(k_Prep_Shims_mTxc.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
end

% FFT back to image space
im_Relativec = fftshift(fftshift(ifft2(ifftshift(ifftshift(k_Relativec,1),2)),1),2);
im_Ref_Shims_mTxc = fftshift(fftshift(ifft2(ifftshift(ifftshift(k_Ref_Shims_mTxc,1),2)),1),2);
im_Prep_Shims_mTxc = fftshift(fftshift(ifft2(ifftshift(ifftshift(k_Prep_Shims_mTxc,1),2)),1),2);

% Combine sensitivity maps with images, sum across receive
% channels without acceleration
Relative_Images = squeeze(sum(bsxfun(@times,im_Relativec,conj(rx_sens(:,:,:,1))),3));
Ref_Shims = squeeze(sum(bsxfun(@times,im_Ref_Shims_mTxc,conj(rx_sens(:,:,:,1))),3));
Prep_Shims = squeeze(sum(bsxfun(@times,im_Prep_Shims_mTxc,conj(rx_sens(:,:,:,1))),3));

plotlowressyntheticimages(Relative_Images,Ref_Shims,Prep_Shims) % Supporting Figure 1

for iter_n  = 1:size(niters,2)
    niters = niters(iter_n); clearvars -except accelerations imaging_plane sz sx sy calib niters iter_n niters_array masks NRepeats Abs_Maps_Recon use_extra_Tx_Voltages image_mask acceleration Enc_Mat Shim_Setting1 Shim_Setting2 location filename Folder1 Folder2 k_Relativec k_Ref_Shims_mTxc k_Prep_Shims_mTxc im_Relativec im_Ref_Shims_mTxc im_Prep_Shims_mTxc Relative_Images Ref_Shims Prep_Shims rx_sens tx_sens kernel eig_thresh Hanning2D calib1 calib2
    load([Folder1,filename,'.mat']); % load in reconstructed kspace data
    filename2 = [filename,'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat']; % filename for reconstructed image data
    if exist([Folder2,filename2],'file') == 2
        disp('Data already exists. Loading in')
        load([Folder2,filename2]);
    else
        % pre-allocate
        err_image_rel = zeros([NRepeats,3,size(accelerations,2)]);
        err_image_ref = zeros([NRepeats,3,size(accelerations,2)]);
        err_image_prep = zeros([NRepeats,3,size(accelerations,2)]);
        err_RelB1 = zeros([NRepeats,3,size(accelerations,2)]);
        err_AbsB1 = zeros([NRepeats,3,size(accelerations,2)]);
        err_FinalB1 = zeros([NRepeats,3,size(accelerations,2)]);
        Rel_Maps = zeros([sz, 8]);
        Rel_Maps_acc = zeros([sz, 8]);
        Maps = zeros([sz, 8, NRepeats,3,size(accelerations,2)]);
        Maps_acc = zeros([sz, 8, NRepeats,3,size(accelerations,2)]);
        rx_sens_acc = zeros([sz, 8, 1, NRepeats,3,size(accelerations,2)]);
        tx_sens_acc = zeros([sz, 8, 1, NRepeats,3,size(accelerations,2)]);
        
        clear B1_Encoded_ref B1_Encoded_acc B1_UnEncoded_ref B1_UnEncoded_acc
        parfor accel_ind = 1:size(accelerations,2)
            for same_masks = 2%1:3 %Code used for testing different alternative masking strategies 
                for repeat_n = 1:NRepeats
                    % Calculate sensitivity maps from cropped k-space for
                    % accelerated data only from Relative images
                    rx_sens_acc(:,:,:,:,repeat_n,same_masks,accel_ind) = tx_espirit(permute(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind),[1,2,4,3]), sz, kernel, eig_thresh);
                    tx_sens_acc(:,:,:,:,repeat_n,same_masks,accel_ind) = tx_espirit(recon_rel(calib1,calib2,:,:,repeat_n,same_masks,accel_ind), sz, kernel, eig_thresh);
                end
            end
        end
        nrmse(rx_sens_acc(:,:,:,:,1,2,1),rx_sens)
        nrmse(rx_sens_acc(:,:,:,:,1,2,3),rx_sens)
        
        if (sz(1) ~= sx || sz(2) ~= sy) && any(size(recon_rel,[1,2]) ~= sz)
            % Hanning filter to padd k-space
            recon_rel = padarray(padarray(recon_rel.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
            recon_ref = padarray(padarray(recon_ref.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
            recon_prep = padarray(padarray(recon_prep.*Hanning2D,ceil((sz-[sx,sy])/2),'pre'),floor((sz-[sx,sy])/2),'post');
        end
        
                    % FFT back to image space
                    im_Rel_acc = fftshift(fftshift(ifft2(ifftshift(ifftshift(recon_rel,1),2)),1),2);
                    im_Ref_acc = fftshift(fftshift(ifft2(ifftshift(ifftshift(recon_ref,1),2)),1),2);
                    im_Prep_acc = fftshift(fftshift(ifft2(ifftshift(ifftshift(recon_prep,1),2)),1),2);
                    
                    % Combine sensitivity maps with images, sum across receive
                    % channels with acceleration
                    Relative_Images_acc = squeeze(sum(bsxfun(@times,im_Rel_acc,conj(rx_sens_acc)),3));
                    Ref_Shims_acc = squeeze(sum(bsxfun(@times,im_Ref_acc,conj(rx_sens_acc)),3));
                    Prep_Shims_acc = squeeze(sum(bsxfun(@times,im_Prep_acc,conj(rx_sens_acc)),3));
                    
                    % Calculate relative maps from transmit sensitivities
                    % or relative images
                    Rel_Maps = squeeze(tx_sens./sum(abs(tx_sens),3));
                    Rel_Maps_acc = squeeze(tx_sens_acc./sum(abs(tx_sens_acc),3));
                    %Rel_Maps = squeeze(Relative_Images./sum(abs(Relative_Images),3));
                    %Rel_Maps_acc = squeeze(Relative_Images_acc./sum(abs(Relative_Images_acc),3));
                    %figure(); imagesc(imtile(abs(Rel_Maps))); figure(); imagesc(imtile(abs(Rel_Maps_acc))); figure(); imagesc(imtile(abs(Rel_Maps)-abs(Rel_Maps_acc)))
                    
                    % Absolute shim maps (un-accelerated)
                    FA_Shim1 = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims(:,:,1:2:end) ./ Ref_Shims(:,:,1:2:end)))));
                    FA_Shim2 = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims(:,:,2:2:end) ./ Ref_Shims(:,:,2:2:end)))));
                    % Absolute shim maps (accelerated)
                    FA_Shim1_acc = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims_acc(:,:,1:2:end,:,:,:) ./ Ref_Shims_acc(:,:,1:2:end,:,:,:)))));
                    FA_Shim2_acc = squeeze((180/pi).*image_mask.*(acos(real(Prep_Shims_acc(:,:,2:2:end,:,:,:) ./ Ref_Shims_acc(:,:,2:2:end,:,:,:)))));
                    %figure(); imagesc(imtile({abs(FA_Shim1),abs(FA_Shim2)}))
                    
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
                    C_acc = ((c1_acc.*abs(Rel_Shim1_acc).^m) + (c2_acc.*abs(Rel_Shim2_acc).^m))./ (abs(Rel_Shim1_acc).^m + abs(Rel_Shim2_acc).^m);
                    C(isnan(C)) = 0; C_acc(isnan(C_acc)) = 0; % Remove NaNs
                    
                    % Scale relative maps to absolute maps
                    Maps = image_mask.*C.*Rel_Maps;
                    Maps_acc = image_mask.*repmat(permute(C_acc,[1 2 6 3 4 5]),[1 1 8]).*Rel_Maps_acc;

        % Calculate difference in absolute B1 maps accelerated to non-accelerated
        save([Folder2,filename2],'Maps','Maps_acc');
    end
    disp(['Iteration ',num2str(niters(iter_n)), ' completed for all acceleration factors ', datestr(now, 'dd/mm/yy-HH:MM')])
end
toc
disp('Finished Simulations.')
end