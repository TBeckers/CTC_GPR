% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function [out] = rownorm(matrix)
out=sqrt(sum(matrix.^2,2));

end

