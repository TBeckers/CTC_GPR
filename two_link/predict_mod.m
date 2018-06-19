% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function [out] = predict_mod(mod,hyp,input_data,test_data)
if mod == 1
    ks=kernel_se(test_data,input_data, hyp.sd,hyp.ell);
    out=kernel_se(test_data,input_data, hyp.sd,hyp.ell)'*(hyp.Ky)+diag(kernel_se(test_data,test_data, hyp.sd,hyp.ell)-ks'*hyp.K*ks)*randn;
elseif mod == 2
    out=kernel_se(test_data,input_data, hyp.sd,hyp.ell)'*(hyp.Ky);
elseif mod == 3
    ks=kernel_se(test_data,input_data, hyp.sd,hyp.ell);
    out=diag(kernel_se(test_data,test_data, hyp.sd,hyp.ell)-ks'*hyp.K*ks);
end
        
end

