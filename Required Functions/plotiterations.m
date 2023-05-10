function [outputArg1,outputArg2] = plotiterations(Maps_acc,niters_array,acceleration)
ROI_xcords = 11:24; ROI_ycords = 12:26;

for iter_n  = 1:size(niters_array,2)
niters = niters_array(iter_n);

for accel_ind = 2:size(acceleration,2)
    for n = 1:3
iter_data(:,n,accel_ind,iter_n) = nonzeros(abs(Maps_acc(ROI_xcords,ROI_ycords,:,:,n,accel_ind) - Maps_acc(ROI_xcords,ROI_ycords,:,:,n,1))); % Crop to ROI and remove nonzeros
rmseiter_data(:,n,accel_ind,iter_n) = nonzeros((Maps_acc(ROI_xcords,ROI_ycords,:,:,n,accel_ind) - Maps_acc(ROI_xcords,ROI_ycords,:,:,n,1))); % Crop to ROI and remove nonzeros
    end
end
end
rmseiter_data = squeeze(sqrt(mean(rmseiter_data.*conj(rmseiter_data),1)));
plotColors = jet(size(acceleration,2));

fig = figure('color','w','Units','centimeters','Position',[0,0,10.526597938144329,8.876907216494844]); tiledlayout('flow','Padding','none','TileSpacing','compact'); 
for accel_n = 1:size(acceleration,2)
plot(niters_array,squeeze(rmseiter_data(2,accel_n,:)),'color',plotColors(accel_n,:)); hold on
end
lgd = legend(split(num2str(acceleration)),'NumColumns',3); lgd.Location = 'northeast';
lgd.Title.String = 'Acc. Factor';
xlabel('Iteration #','Fontsize',axis_fontsize)
ylabel(['RMSE in B_1^+ [',char(176),']'],'Fontsize',axis_fontsize)
xlim([1 niters_array(end)])
ylim([0 20])

end

