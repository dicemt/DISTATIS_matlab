function distData = distatis_s2d(data)
% Copyright (C) 2020 Daisuke MATSUYOSHI

[nrow,ncol] = size(data);
distData = zeros(nrow,nrow,ncol);
for i=1:ncol
    distData(:,:,i) = sort2distance(data(:,i));
end

function distMat = sort2distance(data)
nD = numel(data);
expndata = repmat(data,1,nD);
distMat = 1 - (expndata == transpose(expndata));
