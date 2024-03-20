function [Scaled_Mask] = ScaleMask(Mask,New_Size)
% Reduce size of image mask by FFT and cropping in k-space
k_mask = fft2c(Mask);
min_ind_kz1 = round(round((size(k_mask,1)+1)/2) - New_Size(1)/2);
min_ind_kz2 = round(round((size(k_mask,2)+1)/2) - New_Size(2)/2);
k_mask = k_mask(min_ind_kz1:min_ind_kz1+New_Size(1)-1,min_ind_kz2:min_ind_kz2+New_Size(2)-1);
Scaled_Mask = ifft2c(k_mask);
% Threshold to binarise re-scaled image_mask
Scaled_Mask(abs(Scaled_Mask) < prctile(abs(Scaled_Mask),(1-sum(Mask,'all')./(size(Mask,1)*size(Mask,2)))*1e2,'all')) = 0;
Scaled_Mask(abs(Scaled_Mask) >= prctile(abs(Scaled_Mask),(1-sum(Mask,'all')./(size(Mask,1)*size(Mask,2)))*1e2,'all')) = 1;
end

