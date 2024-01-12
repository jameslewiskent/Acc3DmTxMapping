function plotfigure2(settings)
cmin = 0; cmax = 150;
err_min = -10; err_max = 10;
same_masks = 2; % All Tx, ref/prep images use different masks
axis_fontsize = 11;
ROI_xcords = 14:26; ROI_ycords = 12:26; % ROIs for sz = [32 32]

Recon_Type = settings.Recon_Type;
sx = settings.sx;
sy = settings.sy;
sz = settings.sz;
calib = settings.calib;
niters = settings.niters_array;
accelerations = settings.accelerations;
NRepeats = settings.NRepeats;
Shim_Setting1 = settings.Shim_Setting1;
type = settings.type;

filename = ['Simulated_',Recon_Type,'Recon_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats),'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat'];
load(['Data',filesep,'Synthetic Body Simulation Results',filesep,'ReconData',filesep,filename],'Maps','Maps_acc');

% Adjust ROI incase sz has been altered from [32 32]
ROI_xcords = round((ROI_xcords(1)./32)*sz(1)):round((ROI_xcords(end)./32)*sz(1));
ROI_ycords = round((ROI_ycords(1)./32)*sz(2)):round((ROI_ycords(end)./32)*sz(2));

if strcmp(type,'long')
    figure('color','w','Units','centimeters','Position',[10.710333333333333,-6.731,8.614833333333335,25.5905]); tiledlayout(4,1,'Padding','none','TileSpacing','compact');
else
    figure('color','w','Units','centimeters','Position',[0,0,19.57916666666667,15.425208333333334]); tiledlayout(2,2,'Padding','none','TileSpacing','compact');
end
ax1 = nexttile;
MeanData_to_plot = imtile(reshape(squeeze(mean(abs(permute(Maps_acc(:,:,:,:,same_masks,1:size(accelerations,2)),[1 2 6 4 5 3])),4)),size(Maps_acc,1),size(Maps_acc,2),[]),'GridSize',[size(Maps_acc,3),size(accelerations,2)]);
imagesc(MeanData_to_plot,[cmin cmax]);
cb = colorbar;
cb.Label.String = ['Mean \alpha [',char(176),']']; cb.FontSize = axis_fontsize-2; cb.Label.FontSize = axis_fontsize;
colormap(ax1,turbo)
xticks(((1:size(accelerations,2)).*sz(2)) - sz(2)/2);
xticklabels(strsplit(num2str(accelerations)));
xlabel('Undersampling Factor','FontSize',axis_fontsize)
ylabel('Transmit Channel','FontSize',axis_fontsize)
yticks(((1:size(Maps_acc,3)).*sz(1)) - sz(1)/2);
yticklabels(1:size(Maps_acc,3));
title('a) Mean','Fontsize',axis_fontsize+1);
hold on
rectangle('Position',[ROI_ycords(1)-1,ROI_xcords(1)-1,ROI_ycords(end)-ROI_ycords(1)+2,ROI_xcords(end)-ROI_xcords(1)+2],'LineWidth',0.5,'LineStyle','-','EdgeColor','r')

ax2 = nexttile; % SD
SDData_to_plot = imtile(reshape(squeeze(std(abs(permute(Maps_acc(:,:,:,:,same_masks,1:size(accelerations,2)),[1 2 6 4 5 3])),[],4)),size(Maps_acc,1),size(Maps_acc,2),[]),'GridSize',[size(Maps_acc,3),size(accelerations,2)]);
imagesc(SDData_to_plot,[0 10]);
cb = colorbar;
cb.Label.String = ['SD in \alpha [',char(176),']']; cb.FontSize = axis_fontsize-2; cb.Label.FontSize = axis_fontsize;
colormap(ax2,bluewhitered)
xticks(((1:size(accelerations,2)).*sz(2)) - sz(2)/2);
xticklabels(strsplit(num2str(accelerations)));
xlabel('Undersampling Factor','FontSize',axis_fontsize)
ylabel('Transmit Channel','FontSize',axis_fontsize)
yticks(((1:size(Maps_acc,3)).*sz(1)) - sz(1)/2);
yticklabels(1:size(Maps_acc,3));
title('b) SD','Fontsize',axis_fontsize+1);

% Plot difference in Maps in ROI
Maps_acc_crop = Maps_acc(ROI_xcords,ROI_ycords,:,:,:,1:size(accelerations,2)); Maps_crop = Maps(ROI_xcords,ROI_ycords,:,:,:,:);
MeanData_to_plot = imtile(reshape(squeeze(mean(permute(abs(Maps_acc_crop(:,:,:,:,same_masks,1:size(accelerations,2))) - abs(Maps_crop),[1 2 6 4 5 3]),4)),size(Maps_acc_crop,1),size(Maps_acc_crop,2),[]),'GridSize',[size(Maps_acc_crop,3),size(Maps_acc_crop,6)]);

Maps_crop_CP = sum(bsxfun(@times,Maps_crop,permute(Shim_Setting1,[1 3 2])),3);
Maps_acc_crop_CP = sum(bsxfun(@times,Maps_acc_crop,permute(Shim_Setting1,[1 3 2])),3);
% figure();imagesc(abs(Maps_acc_crop_CP(:,:,1,1,2,1)))
% figure();imagesc(abs(Maps_crop_CP(:,:,1,1,1,1)))

Mean_Maps_crop_CP = mean(abs(nonzeros(Maps_crop_CP)),[1:4]);
SD_Maps_crop_CP = std(abs(nonzeros(Maps_crop_CP)),[],[1:4]);
CV_Maps_crop_CP = squeeze(SD_Maps_crop_CP./Mean_Maps_crop_CP).*100;
disp(['Mean = ',num2str(Mean_Maps_crop_CP)]);
disp(['Min = ',num2str(min(abs(nonzeros(Maps_crop_CP))))]);
disp(['Max = ',num2str(max(abs(nonzeros(Maps_crop_CP))))]);
disp(['CV = ',num2str(CV_Maps_crop_CP)]);

ax3 = nexttile;
imagesc(MeanData_to_plot,[err_min err_max]);
cb = colorbar;
cb.Label.String = ['Difference in \alpha [',char(176),']']; cb.FontSize = axis_fontsize-2; cb.Label.FontSize = axis_fontsize;
colormap(ax3,bluewhitered)
xticks(((1:size(accelerations,2)).*size(Maps_acc_crop,2)) - size(Maps_acc_crop,2)/2);
xticklabels(strsplit(num2str(accelerations)));
xlabel('Undersampling Factor','Fontsize',axis_fontsize)
ylabel('Transmit Channel','Fontsize',axis_fontsize)
yticks(((1:size(Maps_acc_crop,3)).*size(Maps_acc_crop,1)) - size(Maps_acc_crop,1)/2);
yticklabels(1:size(Maps_acc_crop,3));
title('c) Difference in ROI','Fontsize',axis_fontsize+1);

% plot RMSE in B1 Maps
ylimits = [0 10];
datacp = Maps_acc_crop_CP - Maps_crop_CP;
for accel_ind = 1:size(accelerations,2)
    for n = 2
        for repeat_n = 1:size(datacp,4)
            data2cp(:,repeat_n,n,accel_ind) = nonzeros(datacp(:,:,:,repeat_n,n,accel_ind));
            rmsedatacp(repeat_n,n,accel_ind) = squeeze(sqrt(mean(data2cp(:,repeat_n,n,accel_ind).*conj(data2cp(:,repeat_n,n,accel_ind)),1)));
        end
    end
end
nexttile;
warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout') % boxplot bug?
%Markersize = 7;
boxplotgroup = rmsedatacp(:,n,1);
for accel_ind = 2:size(accelerations,2)
    %     %plot(accelerations(accel_ind),mean(abs(data2cp(:,:,2,accel_ind)),1),'k.','MarkerSize',Markersize); hold on
    %     plot(accelerations(accel_ind),(rmsedatacp(:,:,2,accel_ind)),'b.','MarkerSize',Markersize); hold on
    boxplotgroup = cat(2,boxplotgroup,rmsedatacp(:,2,accel_ind));
end
boxplot(boxplotgroup,'Positions',accelerations); hold on
%plot(accelerations,squeeze(mean(rmsedatacp(:,n,:),1)),'x')
xticks(accelerations)
xlim([0 accelerations(end)+1])
ylim(ylimits)
xticklabels(strsplit(num2str(accelerations)))
xlabel('Undersampling Factor','Fontsize',axis_fontsize)
ylabel(['Error in B_1^+ [',char(176),']'],'Fontsize',axis_fontsize)
title('d) CP^+ Error in ROI','Fontsize',axis_fontsize+1);
%lgd = legend('MAE','RMSE'); lgd.Location = 'northwest';
grid on

end

