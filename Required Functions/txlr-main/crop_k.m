function k = crop_k(k, dims)

N1  = ceil((size(k,1)+1)/2);
k   = k(N1-dims(1)/2:N1+dims(2)/2-1,:,:,:);

N2  = ceil((size(k,2)+1)/2);
k   = k(:,N2-dims(2)/2:N2+dims(2)/2-1,:,:);