function plotsyntheticheartmaps(Maps,Maps_acc,sz,acceleration,Shim_Setting1)
cmin = 0; cmax = 130;
err_min = -10; err_max = 10;
same_masks = 2; % All Tx, ref/prep images use different masks
axis_fontsize = 11;
ROI_xcords = 11:24; ROI_ycords = 12:26;


fig = figure('color','w','Units','centimeters','Position',[0,0,19.57916666666667,15.425208333333334]); tiledlayout('flow','Padding','none','TileSpacing','compact'); 
ax1 = nexttile;
MeanData_to_plot = imtile(reshape(squeeze(mean(abs(permute(Maps_acc(:,:,:,:,same_masks,:),[1 2 6 4 5 3])),4)),size(Maps_acc,1),size(Maps_acc,2),[]),'GridSize',[size(Maps_acc,3),size(Maps_acc,6)]);
imagesc(MeanData_to_plot,[0 180]);
cb = colorbar;
cb.Label.String = ['Mean \alpha [',char(176),']']; cb.FontSize = axis_fontsize-2; cb.Label.FontSize = axis_fontsize;
colormap(ax1,turbo)
xticks(((1:size(acceleration,2)).*sz(2)) - sz(2)/2);
xticklabels(strsplit(num2str(acceleration)));
xlabel('Undersampling Factor','FontSize',axis_fontsize)
ylabel('Transmit Channel','FontSize',axis_fontsize)
yticks(((1:size(Maps_acc,3)).*sz(1)) - sz(1)/2);
yticklabels(1:size(Maps_acc,3));
title(['a) Mean'],'Fontsize',axis_fontsize+1);
hold on
rectangle('Position',[ROI_ycords(1)-1,ROI_xcords(1)-1,ROI_ycords(end)-ROI_ycords(1)+2,ROI_xcords(end)-ROI_xcords(1)+2],'LineWidth',0.5,'LineStyle','-','EdgeColor','r')

ax2 = nexttile; % SD
SDData_to_plot = imtile(reshape(squeeze(std(abs(permute(Maps_acc(:,:,:,:,same_masks,:),[1 2 6 4 5 3])),[],4)),size(Maps_acc,1),size(Maps_acc,2),[]),'GridSize',[size(Maps_acc,3),size(Maps_acc,6)]);
imagesc(SDData_to_plot,[0 10]);
cb = colorbar;
cb.Label.String = ['SD in \alpha [',char(176),']']; cb.FontSize = axis_fontsize-2; cb.Label.FontSize = axis_fontsize;
colormap(ax2,bluewhitered)
xticks(((1:size(acceleration,2)).*sz(2)) - sz(2)/2);
xticklabels(strsplit(num2str(acceleration)));
xlabel('Undersampling Factor','FontSize',axis_fontsize)
ylabel('Transmit Channel','FontSize',axis_fontsize)
yticks(((1:size(Maps_acc,3)).*sz(1)) - sz(1)/2);
yticklabels(1:size(Maps_acc,3));
title(['b) SD'],'Fontsize',axis_fontsize+1);

% Plot difference in Maps in ROI
Maps_acc_crop = Maps_acc(ROI_xcords,ROI_ycords,:,:,:,:); Maps_crop = Maps(ROI_xcords,ROI_ycords,:,:,:,:); 
MeanData_to_plot = imtile(reshape(squeeze(mean(permute(abs(Maps_acc_crop(:,:,:,:,same_masks,:)) - abs(Maps_crop),[1 2 6 4 5 3]),4)),size(Maps_acc_crop,1),size(Maps_acc_crop,2),[]),'GridSize',[size(Maps_acc_crop,3),size(Maps_acc_crop,6)]);

Maps_crop_CP = sum(bsxfun(@times,Maps_crop,permute(Shim_Setting1,[1 3 2])),3);
Maps_acc_crop_CP = sum(bsxfun(@times,Maps_acc_crop,permute(Shim_Setting1,[1 3 2])),3);

Mean_Maps_crop_CP = mean(abs(nonzeros(Maps_crop_CP)),[1:4])
SD_Maps_crop_CP = std(abs(nonzeros(Maps_crop_CP)),[],[1:4]);
CV_Maps_crop_CP = squeeze(SD_Maps_crop_CP./Mean_Maps_crop_CP).*100
min(abs(nonzeros(Maps_crop_CP)))
max(abs(nonzeros(Maps_crop_CP)))

ax3 = nexttile;
imagesc(MeanData_to_plot,[err_min err_max]);
cb = colorbar;
cb.Label.String = ['Difference in \alpha [',char(176),']']; cb.FontSize = axis_fontsize-2; cb.Label.FontSize = axis_fontsize;
colormap(ax3,bluewhitered)
xticks(((1:size(acceleration,2)).*size(Maps_acc_crop,2)) - size(Maps_acc_crop,2)/2);
xticklabels(strsplit(num2str(acceleration)));
xlabel('Undersampling Factor','Fontsize',axis_fontsize)
ylabel('Transmit Channel','Fontsize',axis_fontsize)
yticks(((1:size(Maps_acc_crop,3)).*size(Maps_acc_crop,1)) - size(Maps_acc_crop,1)/2);
yticklabels(1:size(Maps_acc_crop,3));
title(['c) Difference in ROI'],'Fontsize',axis_fontsize+1);

% plot MAD in B1 Maps
Markersize = 7;
ylimits = [0 14];

data = (Maps_acc_crop) - (Maps_acc_crop(:,:,:,:,:,1)); 
datacp = (Maps_acc_crop_CP) - (Maps_acc_crop_CP(:,:,:,:,:,1)); 

for accel_ind = 2:size(acceleration,2)
    for n = 2
data2cp(:,n,accel_ind) = nonzeros(datacp(:,:,:,:,n,accel_ind));
rmsedatacp(:,n,accel_ind) = sqrt(mean((data2cp(:,n,accel_ind)).*conj(data2cp(:,n,accel_ind)),1)); 
    end
end

ax4 = nexttile;
for accel_ind = 1:size(acceleration,2)
    plot(acceleration(accel_ind),mean(abs(data2cp(:,2,accel_ind)),1),'k.','MarkerSize',Markersize); hold on
    plot(acceleration(accel_ind),(rmsedatacp(:,2,accel_ind)),'b.','MarkerSize',Markersize); hold on
end
xticks(acceleration)
xlim([1 acceleration(end)+1])
ylim(ylimits)
xticklabels(strsplit(num2str(acceleration)))
xlabel('Undersampling Factor','Fontsize',axis_fontsize)
ylabel(['Error in B_1^+ [',char(176),']'],'Fontsize',axis_fontsize)
title(['d) CP^+ Error in ROI'],'Fontsize',axis_fontsize+1);

lgd = legend('MAE','RMSE'); lgd.Location = 'northwest';
grid on

end

