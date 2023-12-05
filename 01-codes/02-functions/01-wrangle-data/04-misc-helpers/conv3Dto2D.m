% This function converts a 3D matrix to 2D ensuring is concatenates by row.
% Otherwise simply using reshape will concatenate along columns.
function out = conv3Dto2D(in)

out = reshape(permute(in,[1,3,2]),[],size(in,2),1);

end