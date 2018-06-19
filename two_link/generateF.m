% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function [ out] = generateF(f1,f2,data_vec)
[q1,q2,dq1,dq2]=ndgrid(data_vec);
q=[q1(:),q2(:)];
dq=[dq1(:),dq2(:)];

training_input=[dq,q];
training_output1=f1(dq,q);
training_output2=f2(dq,q);


likfunc = @likGauss;
covfunc = @covSEard;
    
    hyp.cov = ones(5,1); hyp.lik = log(1);
    hyp = minimize(hyp, @gp, -200, @infExact, [], covfunc, likfunc, training_input, training_output1);
    
    sd1=sqrt(exp(2*hyp.cov(end)));
    ell1=exp(hyp.cov(1:end-1));
    
    K1=kernel_se(training_input',training_input', sd1,ell1)+eye(length(training_input))*exp(2*hyp.lik);
    
    hyp.cov = ones(5,1); hyp.lik = log(1);
    hyp = minimize(hyp, @gp, -200, @infExact, [], covfunc, likfunc, training_input, training_output2);
    
    sd2=sqrt(exp(2*hyp.cov(end)));
    ell2=exp(hyp.cov(1:end-1));
    
    K2=kernel_se(training_input',training_input', sd2,ell2)+eye(length(training_input))*exp(2*hyp.lik);
    temp1=K1\training_output1;
    temp2=K2\training_output2;
    out=@(dq,q) 2*[kernel_se([dq;q],training_input', sd1,ell1)'*temp1 ;kernel_se([dq;q],training_input', sd2,ell2)'*temp2];
end

