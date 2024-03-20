function [c] = rmse(a,b,type)
if nargin <3
    type = 'notnorm';
end


if strcmp(type,'norm')
        c   =  norm(a(:)-b(:))/norm(b(:));
    else
        c = sqrt(mean((a - b).*conj(a - b),'all'));
    end
end
