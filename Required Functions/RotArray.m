function [TwoDArray] = RotArray(OneDArray)
N = length(OneDArray);
M = (N-1) / 2;
n = 2 / M * (-M:M);
[x,y] = meshgrid(n);
r = sqrt( x.^2 + y.^2 );
TwoDArray = zeros(N);
TwoDArray(:) = interp1(n, OneDArray, r(:));
TwoDArray(isnan(TwoDArray)) = 0;
end

