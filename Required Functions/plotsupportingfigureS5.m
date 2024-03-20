function plotsupportingfigureS5(settings)
ROI_xcords = 14:26; ROI_ycords = 12:26; % ROIs for sz = [32 32]
axis_fontsize = 11;

Recon_Type = settings.Recon_Type;
sx = settings.sx;
sy = settings.sy;
sz = settings.sz;
calib = settings.calib;
NRepeats = settings.NRepeats;
niters_array = settings.niters_array;
accelerations = settings.accelerations;
Shim_Setting1 = settings.Shim_Setting1;

% Adjust incase sz has been altered from [32 32]
ROI_xcords = round((ROI_xcords(1)./32)*sz(1)):round((ROI_xcords(end)./32)*sz(1));
ROI_ycords = round((ROI_ycords(1)./32)*sz(2)):round((ROI_ycords(end)./32)*sz(2));

for iter_n  = 1:size(niters_array,2)
    niters = niters_array(iter_n);
    filename = ['Simulated_',Recon_Type,'Recon_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats),'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat'];
    load(['Data',filesep,'Synthetic Body Simulation Results',filesep,'ReconData',filesep,filename],'Maps','Maps_acc');
    
    Maps_CP = sum(bsxfun(@times,Maps,permute(Shim_Setting1,[1 3 2])),3);
    Maps_acc_CP = sum(bsxfun(@times,Maps_acc,permute(Shim_Setting1,[1 3 2])),3);
    
    for accel_ind = 1:size(accelerations,2)
        for n = 1:3
            for repeat_n = 1:NRepeats
                diffiter_data(:,repeat_n,n,accel_ind,iter_n) = nonzeros((Maps_acc_CP(ROI_xcords,ROI_ycords,:,repeat_n,n,accel_ind) - Maps_CP(ROI_xcords,ROI_ycords))); % Crop to ROI and remove nonzeros
            end
        end
    end
end
rmseiter_data = squeeze(sqrt(mean(diffiter_data.*conj(diffiter_data),1)));


plotColors = jet(size(accelerations,2));

figure('color','w','Units','centimeters','Position',[13.059833333333334,2.963333333333333,13.631333333333334,11.557]); tiledlayout('flow','Padding','none','TileSpacing','compact');
for accel_n = 1:size(accelerations,2)
    plot(niters_array,squeeze(mean(rmseiter_data(:,2,accel_n,:),1)),'color',plotColors(accel_n,:),'LineWidth',1.5); hold on
    plot(niters_array,squeeze(mean(rmseiter_data(:,2,accel_n,:),1))+1.96*squeeze(std(rmseiter_data(:,2,accel_n,:),[],1)),'color',plotColors(accel_n,:),'LineStyle','--','HandleVisibility','off');
    plot(niters_array,squeeze(mean(rmseiter_data(:,2,accel_n,:),1))-1.96*squeeze(std(rmseiter_data(:,2,accel_n,:),[],1)),'color',plotColors(accel_n,:),'LineStyle','--','HandleVisibility','off');
    
end
lgd = legend(split(num2str(accelerations)),'NumColumns',3); lgd.Location = 'northeast';
lgd.Title.String = 'Acc. Factor';
xlabel('Iteration #','Fontsize',axis_fontsize)
ylabel(['CP^+ Error in ROI [',char(176),']'],'Fontsize',axis_fontsize)
xlim([niters_array(1) niters_array(end)])
ylim([0 20])
end

