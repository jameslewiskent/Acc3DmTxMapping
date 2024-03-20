function plotsupportingfigureS1(Relative_Images,Ref_Shims,Prep_Shims)
% Just for plotting synthetic low resolution images (without worrying about receive channels for figure)
count = 1;
for n = 1:size(Ref_Shims,3)
catdata(:,:,count) = abs(Ref_Shims(:,:,n)); count = count +1;
catdata(:,:,count) = abs(Prep_Shims(:,:,n)); count = count +1;
end

plotdata = cat(3,abs(Relative_Images),catdata);
figure('color','w','units','centimeters','Position',[10.625666666666667,5.7785,22.37316666666667,12.170833333333334]); tiledlayout('flow','TileSpacing','none','Padding','none'); nexttile;
imagesc(imtile(abs(plotdata)),[0 prctile(abs(plotdata),99,'all')]); axis image off
cb = colorbar;
cb.Label.String = 'Image Magnitude, [a.u.]';
end