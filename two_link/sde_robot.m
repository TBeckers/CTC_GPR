% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

function dxdt = sde_robot(t,x,xp,n,tau,H,C,G,F)

    dxdt=zeros(2*n,1);
    dxdt(1:n)=x(n+1:2*n)-xp(1:n);
    dxdt(n+1:2*n)=H(x(1:n))*xp(n+1:end)+C(x(n+1:2*n),x(1:n))*x(n+1:2*n)+G(x(1:n))+F(xp(n+1:end),x(n+1:2*n),x(1:n))-tau(t,xp(n+1:end),x(n+1:2*n),x(1:n));
end

