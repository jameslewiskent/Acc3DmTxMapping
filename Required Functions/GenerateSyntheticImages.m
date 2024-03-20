function GenerateSyntheticImages(Enc_Scheme)
load(['Data',filesep,'SyntheticDuke.mat'],'B1Tx','B1Rx','PDimage','sigma');

cor_sl = 1:139; sag_sl = 1:178; ax_sl = 60; 
RF_Voltage = 60; Gamma = 42.57747892e6;

% Centre body in FOV
B1Tx = circshift(B1Tx,15,1);
B1Rx = circshift(B1Rx,15,1);
PDimage = circshift(PDimage,15,1);
sigma = circshift(sigma,15,1);

B1Tx = B1Tx(cor_sl,sag_sl,ax_sl,:).*RF_Voltage/sqrt(50);

% Calc Encoding matrix
[Enc_Mat,Modes] = CalcEncMat(Enc_Scheme);

% Calculate B1 Tx
B1Tx_modes = zeros([size(B1Tx,1:2),1,Modes]);
for mode = 1:Modes
    B1Tx_modes(:,:,1,mode) = sum(bsxfun(@times,B1Tx,permute(Enc_Mat(mode,:),[1,3,4,2])),4); % Sum of channels for Slice 60 in (T)
end
Tx_FA_modes = B1Tx_modes.*Gamma.*1e-3.*360; % FA (degrees) (this is complex)

Relative_Rx = bsxfun(@times,bsxfun(@times,PDimage(cor_sl,sag_sl,ax_sl),permute(B1Rx(cor_sl,sag_sl,ax_sl,:),[1,2,4,3])),Tx_FA_modes(:,:,1,1:8));
Ref_Shim1_Rx = bsxfun(@times,bsxfun(@times,PDimage(cor_sl,sag_sl,ax_sl),permute(B1Rx(cor_sl,sag_sl,ax_sl,:),[1,2,4,3])),Tx_FA_modes(:,:,1,9));
Ref_Shim2_Rx = bsxfun(@times,bsxfun(@times,PDimage(cor_sl,sag_sl,ax_sl),permute(B1Rx(cor_sl,sag_sl,ax_sl,:),[1,2,4,3])),Tx_FA_modes(:,:,1,10));

mz_prep = cos(abs(Tx_FA_modes(:,:,1,9:10)*pi/180));
Prep_Shim1_Rx = bsxfun(@times,Ref_Shim1_Rx,mz_prep(cor_sl,sag_sl,1,1));
Prep_Shim2_Rx = bsxfun(@times,Ref_Shim2_Rx,mz_prep(cor_sl,sag_sl,1,2));

% Remove lungs from mask
image_mask_full = ones(size(PDimage(cor_sl,sag_sl,ax_sl)));
image_mask_full(PDimage(cor_sl,sag_sl,ax_sl) <= 80) = 0; % Full size image_mask

% Define Heart Mask
slice = sigma(cor_sl,sag_sl,ax_sl); % For now only consider one slice of conductivity data to mask
mask_heart = zeros(size(slice)); % pre-allocate heart mask
for indi = 1:size(slice,1)
    for indj = 1:size(slice,2)
        if (slice(indi,indj) > 0.85 && slice(indi,indj) < 0.95) || (slice(indi,indj) > 1.3 && slice(indi,indj) < 1.4)
            mask_heart(indi,indj) = 1;
        end
    end
end

% Dilate/erode to remove voids and single voxels
se1 = strel('square',3); se2 = strel('square',5);
mask_heart = imerode(imdilate(mask_heart, se1),se2);
mask_heart = imdilate(mask_heart, se2);
image_heartmask_full = image_mask_full.*mask_heart;


% Concatenate shims together
Ref_Shims = cat(4,Ref_Shim1_Rx,Ref_Shim2_Rx);
Prep_Shims = cat(4,Prep_Shim1_Rx,Prep_Shim2_Rx);

save(['Data',filesep,'SyntheticBodyImages'],'Relative_Rx','Ref_Shims','Prep_Shims','image_mask_full','image_heartmask_full');
end