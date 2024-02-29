function GenerateSyntheticImages(Enc_Scheme)
load(['Data',filesep,'SyntheticDuke.mat'],'B1Tx','B1Rx','PDimage');

cor_sl = 1:139; sag_sl = 1:178; ax_sl = 60;
B1Tx = B1Tx(cor_sl,sag_sl,ax_sl,:);

% Calc Encoding matrix
[Enc_Mat,Modes] = CalcEncMat(Enc_Scheme);

% Calculate B1 Tx
B1Tx_modes = zeros([size(B1Tx,1:2),1,Modes]);
for mode = 1:Modes
    B1Tx_modes(:,:,1,mode) = sum(bsxfun(@times,B1Tx,permute(Enc_Mat(mode,:),[1,3,4,2])),4); % Sum of channels for Slice 60 in (T)
end

mz_prep = cos(abs(B1Tx_modes(:,:,1,9:10)*pi*1e6));

Relative_Rx = bsxfun(@times,bsxfun(@times,PDimage(cor_sl,sag_sl,ax_sl),permute(B1Rx(cor_sl,sag_sl,ax_sl,:),[1,2,4,3])),B1Tx_modes(:,:,1,1:8));
Ref_Shim1_Rx = bsxfun(@times,bsxfun(@times,PDimage(cor_sl,sag_sl,ax_sl),permute(B1Rx(cor_sl,sag_sl,ax_sl,:),[1,2,4,3])),B1Tx_modes(:,:,1,9));
Ref_Shim2_Rx = bsxfun(@times,bsxfun(@times,PDimage(cor_sl,sag_sl,ax_sl),permute(B1Rx(cor_sl,sag_sl,ax_sl,:),[1,2,4,3])),B1Tx_modes(:,:,1,10));
Prep_Shim1_Rx = bsxfun(@times,Ref_Shim1_Rx,mz_prep(cor_sl,sag_sl,1,1));
Prep_Shim2_Rx = bsxfun(@times,Ref_Shim2_Rx,mz_prep(cor_sl,sag_sl,1,2));

% Remove lungs from mask
image_mask_full = ones(size(PDimage(cor_sl,sag_sl,ax_sl)));
image_mask_full(PDimage(cor_sl,sag_sl,ax_sl) <= 80) = 0; % Full size image_mask

Ref_Shims = cat(4,Ref_Shim1_Rx,Ref_Shim2_Rx);
Prep_Shims = cat(4,Prep_Shim1_Rx,Prep_Shim2_Rx);

save(['Data',filesep,'SyntheticBodyImages'],'Relative_Rx','Ref_Shims','Prep_Shims','image_mask_full');
end