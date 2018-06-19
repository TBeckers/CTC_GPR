% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function [ ddq ] = dp(q2,dq,q)
m1=1;
m2=10;
l1=1;
l2=1;
g=1;

ddq=(m2*l1*dq.^2.*sin(q-q2)-m2*g*sin(q2))./(m2*l1*cos(q-q2)-(m1+m2)*l2./cos(q-q2));

end

