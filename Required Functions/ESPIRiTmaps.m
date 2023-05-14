function [maps,weights,cgESPIRiTrecon] = ESPIRiTmaps(calib,kernel,eig_thresh,eig_thresh2,plotyn,DATA)
if nargin < 6
   performCGESPIRiT = 0;
   cgESPIRiTrecon = 'NaN';
else
    performCGESPIRiT = 1;
end

calib = squeeze(calib);

if plotyn == 1
figure(); imagesc(imtile(abs(squeeze(calib))))
end
im = ifft2c(calib);

[sx,sy,Nc] = size(calib);

% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
[k,S] = dat2Kernel(calib,kernel);
idx = max(find(S >= S(1)*eig_thresh));

if plotyn == 1
kdisp = reshape(k,[kernel(1)*kernel(2)*Nc,kernel(1)*kernel(2)*Nc]);
figure, subplot(211), plot([1:kernel(1)*kernel(2)*Nc],S,'LineWidth',2);
hold on;
plot([1:kernel(1)*kernel(2)*Nc],ones(1,kernel(1)*kernel(2)*Nc).*S(1)*eig_thresh,'r-','LineWidth',2);
plot([idx,idx],[0,S(1)],'g--','LineWidth',2)
legend('Singular vector Value','threshold')
title('Singular Vectors')
subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
xlabel('Singular Value #');
title('Singular Vectors')
end

[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

if plotyn == 1
figure, imshow3(abs(W),[],[1,8]);
title('Eigen Values in Image space');
colormap((gray(256))); colorbar;

figure, imshow3(abs(M),[],[8,8]);
title('Magnitude of Eigen Vectors');
colormap(gray(256)); colorbar;

figure, imshow3(angle(M),[],[8,8]);
title('Magnitude of Eigen Vectors');
colormap(jet(256)); colorbar;
end

P = sum(repmat(im,[1,1,1,Nc]).*conj(M),3);

if plotyn == 1
figure, imshow3(abs(P),[],[1,8]);
title('Magnitude of Projection onto Eigen Vectors');
colormap(sqrt(gray(256))); colorbar;

figure, imshow3(angle(P),[],[1,8]);
title('Phase of projection onto Eigen Vectors');
colormap((jet(256))); colorbar;
end

maps = M(:,:,:,end-1:end).*repmat(permute(W(:,:,end-1:end)>eig_thresh2,[1,2,4,3]),[1,1,Nc,1]);

if plotyn == 1
figure, imshow3(abs(maps),[],[2,8]);
title('Absolute sensitivity maps');
colormap((gray(256))); colorbar;

figure, imshow3(angle (maps),[],[2,8]);
title('Phase of sensitivity maps');
colormap((jet(256))); colorbar;
end

% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,end-1:end) ;
weights = (weights - eig_thresh2)./(1-eig_thresh2).* (W(:,:,end-1:end) > eig_thresh2);
weights = -cos(pi*weights)/2 + 1/2;

if performCGESPIRiT ~= 0 % calculate CGESPIRiT

if plotyn == 1
    figure(); imagesc(imtile(abs(weights)))
end

% create ESPIRiT operator
ESP = ESPIRiT(maps,weights);
nIterCG = 12;
for DATA_n = 1:size(DATA,4)
[cgESPIRiTrecon(:,:,:,DATA_n),~] = cgESPIRiT(DATA(:,:,:,DATA_n),ESP, nIterCG, 0.01);
end
end

end

