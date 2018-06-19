% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function [dxdt, u] = sde_robot(t,x,tau,H,C,G,F)

    u=tau(t,0,x(1)+noise,x(2)+noise);

    dxdt=zeros(3,1);
    dxdt(2)=x(1);
    dxdt(1)=H(x(2))\(u-C(x(1),x(2))*x(1)-G(x(2))-F(0,x(1),x(2)));
    
    

