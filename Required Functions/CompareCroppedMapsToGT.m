%% Ground-truth to compare against
clearvars
load(['Data',filesep,'SyntheticBodyImages.mat']','Ref_Shims','Prep_Shims','Relative_Rx','image_mask_full','image_heartmask_full')
Images = cat(4,Relative_Rx,Ref_Shims,Prep_Shims);
sx = 24; sy = 24;
calibsize = [sx sy]; % ESPIRiT calibration size
centre = ceil(size(Images,1:2)/2);
calib1 = centre(1)-floor((calibsize(1)-1)/2):centre(1)+ceil((calibsize(1)-1)/2); % Calculate extent of ESPIRiT calibration region
calib2 = centre(2)-floor((calibsize(2)-1)/2):centre(2)+ceil((calibsize(2)-1)/2); % Calculate extent of ESPIRiT calibration region
kdata = fft2c(Images);
clear Images
%% Recon Ground Truth
rx_sens_GT = tx_espirit(permute(kdata(calib1,calib2,:,1:8),[1,2,4,3]), size(kdata,1:2), [5 5], 0.02);
tx_sens_GT = tx_espirit(kdata(:,:,:,1:8), size(kdata,1:2), [5 5], 0.02);
Rel_Maps_GT = image_mask_full.*squeeze(tx_sens_GT./sum(abs(tx_sens_GT),3));

Images_GT = image_mask_full.*ifft2c(kdata);
Images_GT = squeeze(sum(bsxfun(@times,Images_GT,conj(rx_sens_GT)),3)); % JK removed conj

Images_FA_GT = (180./pi).*acos(real(Images_GT(:,:,11:12)./Images_GT(:,:,9:10)));
Images_FA_GT(isnan(Images_FA_GT)) = 0;

%% Recon cropped images
recon_size = [139 178];
image_mask = ScaleMask(image_mask_full,recon_size);
image_heartmask = ScaleMask(image_heartmask_full,recon_size);

% Define windowing filters
TukeyRot = RotArray(tukeywin(sx+2,0.7)); TukeyRot = TukeyRot(2:end-1,2:end-1);

kdata_c = kdata(calib1,calib2,:,:); % Crop k-space
rx_sens = tx_espirit(permute(kdata_c(:,:,:,1:8),[1,2,4,3]), recon_size, [5 5], 0.02);
tx_sens = tx_espirit(kdata_c(:,:,:,1:8), recon_size, [5 5], 0.02);
Rel_Maps_cropped = image_mask.*squeeze(tx_sens./sum(abs(tx_sens),3));

kdata_p = padarray(padarray(kdata_c.*TukeyRot,ceil((recon_size-[sx,sy])/2),'pre'),floor((recon_size-[sx,sy])/2),'post'); % Zero-pad

Images_cropped = image_mask.*ifft2c(kdata_p);
Images_cropped = squeeze(sum(bsxfun(@times,Images_cropped,conj(rx_sens)),3));
%Images_cropped = squeeze(sum(Images_cropped,3));

Images_FA_cropped = (180./pi).*acos(real(Images_cropped(:,:,11:12)./Images_cropped(:,:,9:10)));

Images_FA_cropped(isnan(Images_FA_cropped)) = 0;


%% Plot Absolute Maps

figure('color','w','Position',[600.2,987.4000000000001,680.8,852]); tiledlayout(3,1,'Padding','none','TileSpacing','none')
ax1 = nexttile(); imagesc(imtile(abs(Images_FA_GT),'GridSize',[1,2]),[0 180]); axis image; colormap(ax1,turbo); colorbar; set(gca,'XTick',[], 'YTick', [])
ax1 = nexttile(); imagesc(imtile(abs(Images_FA_cropped),'GridSize',[1,2]),[0 180]); axis image; colormap(ax1,turbo); colorbar; set(gca,'XTick',[], 'YTick', [])

Diff_Data = imtile(abs(Images_FA_cropped)-abs(Images_FA_GT),'GridSize',[1,2]);
Heart_Boundary = imtile(repmat(double(boundarymask(image_heartmask,4)),[1 1 2]),'GridSize',[1,2]);
imAlpha = ones(size(Diff_Data));
imAlpha(Heart_Boundary == 1) = 0;

Rel_Diff_Data = Diff_Data./imtile(repmat(abs(Images_FA_GT),[1 1 2]),'GridSize',[1,2])*100;
Rel_Diff_Data(isnan(Rel_Diff_Data)) = 0;

ax2 = nexttile(); imagesc(Rel_Diff_Data,'AlphaData',imAlpha,[-10 10]); axis image; colormap(ax2,bluewhitered); cb = colorbar;  set(gca,'XTick',[], 'YTick', [])
set(ax2,'color',[0 1 0]);

nrmse(nonzeros((image_heartmask.*Images_FA_cropped)),nonzeros((image_heartmask.*Images_FA_GT)))
rmse(nonzeros((image_heartmask.*Images_FA_cropped)),nonzeros((image_heartmask.*Images_FA_GT)))

%% Plot Relative Maps 

figure('color','w','Position',[131.4,965.8000000000001,1715.2,584.0000000000001]); tiledlayout(3,1,'Padding','none','TileSpacing','none')
ax1 = nexttile(); imagesc(imtile(abs(Rel_Maps_GT),'GridSize',[1,8]),[0 1]); axis image; colormap(ax1,turbo); colorbar; set(gca,'XTick',[], 'YTick', [])
ax1 = nexttile(); imagesc(imtile(abs(Rel_Maps_cropped),'GridSize',[1,8]),[0 1]); axis image; colormap(ax1,turbo); colorbar; set(gca,'XTick',[], 'YTick', [])

Diff_Data = imtile(abs(Rel_Maps_cropped)-abs(Rel_Maps_GT),'GridSize',[1,8]);
Heart_Boundary = imtile(repmat(double(boundarymask(image_heartmask,4)),[1 1 8]),'GridSize',[1,8]);
imAlpha = ones(size(Diff_Data));
imAlpha(Heart_Boundary == 1) = 0;

Rel_Diff_Data = Diff_Data./imtile(repmat(abs(Rel_Maps_GT),[1 1 8]),'GridSize',[1,8]).*100;
Rel_Diff_Data(isnan(Rel_Diff_Data)) = 0;

ax2 = nexttile(); imagesc(Rel_Diff_Data,'AlphaData',imAlpha,[-10 10]); axis image; colormap(ax2,bluewhitered); cb = colorbar;  set(gca,'XTick',[], 'YTick', [])
set(ax2,'color',[0 1 0]);

nrmse(nonzeros((image_heartmask.*Rel_Maps_cropped)),nonzeros((image_heartmask.*Rel_Maps_GT)))
rmse(nonzeros((image_heartmask.*Rel_Maps_cropped)),nonzeros((image_heartmask.*Rel_Maps_GT)))
