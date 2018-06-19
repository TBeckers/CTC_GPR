% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function [ para ] = odepassout( func,t,x,n_para)
n=length(t);
para=zeros(n,n_para);
for i=1:n
    [~,para(i,:)]=func(t(i),x(i,:));
end

end

