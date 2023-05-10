function [nrmse] = nrmse(a,b)
nrmse = norm(a(:)-b(:))/norm(b(:));
end

