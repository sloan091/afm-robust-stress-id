% This function converts a 4D matrix to 2D ensuring is concatenates by row.
% Otherwise simply using reshape will concatenate along columns.
function out = conv4Dto2D(in)
out = conv3Dto2D(squeeze(reshape(permute(in,[1,3,2,4]),[],size(in,2),...
    1,size(in,4))));
end