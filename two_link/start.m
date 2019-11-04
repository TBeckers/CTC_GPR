% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018
%
% Perfect tracking control for real-world Euler-Lagrange systems is challenging due to uncertainties in the system model 
% and external disturbances. The magnitude of the tracking error can be reduced either by increasing the feedback gains or 
% improving the model of the system. The latter is clearly preferable as it allows to maintain good tracking performance at 
% low feedback gains. However, accurate models are often difficult to obtain.
% We address the problem of stable high-performance tracking control for unknown Euler-Lagrange systems. In particular, 
% we employ Gaussian Process regression to obtain a data-driven model that is used for the feed-forward compensation of 
% unknown dynamics of the system. The model fidelity is used to adapt the feedback gains allowing low feedback gains in 
% state space regions of high model confidence. The proposed control law guarantees a globally bounded tracking error with
% a specific probability.
% In the following example, the benefit of the CTC-GPR is shown for a 2-link robotic manipulator.

% Robot model
n=2;

m1=1;
m2=1;
l1=1;
l2=1;
g=10;

alpha=m1*(l1/2)^2+m2*(l1^2+(l2/2)^2);
beta=m2*l1*l2/2;
delta=m2*(l2/2)^2;
H=@(q) [alpha+2*beta*cos(q(2)), delta+beta*cos(q(2));delta+beta*cos(q(2)), delta];
C=@(dq,q) [-beta*sin(q(2))*dq(2), -beta*sin(q(2))*(dq(1)+dq(2)); beta*sin(q(2))*dq(1), 0];
G=@(q) [(m1+m2)*g*l1/2*sin(q(1))+m2*g*l2/2*sin(q(1)+q(2)); m2*g*l2/2*sin(q(1)+q(2))];

% Estimated model
m1h=0.9;
m2h=1.1;
l1h=0.9;
l2h=1.1;
gh=10;

alpha=m1h*(l1h/2)^2+m2h*(l1h^2+(l2h/2)^2);
beta=m2h*l1h*l2h/2;
delta=m2h*(l2h/2)^2;

hatH=@(q) [alpha+2*beta*cos(q(2)), delta+beta*cos(q(2));delta+beta*cos(q(2)), delta];
hatC=@(dq,q) [-beta*sin(q(2))*dq(2), -beta*sin(q(2))*(dq(1)+dq(2)); beta*sin(q(2))*dq(1), 0];
hatG=@(q) [(m1+m2)*g*l1/2*sin(q(1))+m2*g*l2/2*sin(q(1)+q(2)); m2*g*l2/2*sin(q(1)+q(2))];

F=@(ddq,dq,q) [sin(2*dq(2,:))+cos(2*q(1,:))+ddq(1,:);sin(2*dq(2,:))+2*sin(dq(1,:))];

% Desired trajectory
qd=@(t) [sin(t);cos(t)];
dqd=@(t) [cos(t);-sin(t)];
ddqd=@(t) [-sin(t);-cos(t)];

% Control law
Kp= 10*eye(2);
Kd= 10*eye(2);

u=@(t,ddq,dq,q) hatH(q)*ddqd(t)+hatC(dq,q)*dqd(t)+hatG(q)-Kp*(q-qd(t))-Kd*(dq-dqd(t));

opts = odeset('AbsTol',1e-3);

%u=@(t,dq,q) [0;0];
%% Simulation with classical CTC
tspan=[0:0.05:20];

q0=[0;1];
dq0=[1;0];
ddq0=[0;0];

xp0=[dq0;ddq0];
x0=[q0;dq0];
[t,x_ct] = ode15i(@(t,x,xp) sde_robot(t,x,xp,n,u,H,C,G,F), tspan,x0,xp0,opts);

qd_vec=qd(t')';
dqd_vec=dqd(t')';

%% Plot
figure(1);
clf
title('Control with classical CTC');
subplot(2,1,1)
plot(t,qd_vec(:,1),'--',t,dqd_vec(:,1),'--');
hold on
plot(t,x_ct(:,[1 3]));
legend('qd1','dqd1','q1','dq1');
xlabel('time');
ylabel('States');

subplot(2,1,2)
plot(t,qd_vec(:,2),'--',t,dqd_vec(:,2),'--');
hold on
plot(t,x_ct(:,[2 4]));
legend('qd2','dqd2','q2','dq2');
xlabel('time');
ylabel('States');

drawnow();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data generation for learning (direct q,dq,ddq)
disp('Starting training data generation');
data_vec=0:0.3:1;
data_vec1=-1:1:1;
data_vec2=0:1:1;
[q1,q2,dq1,dq2,ddq1,ddq2]=ndgrid(data_vec,data_vec,data_vec1,data_vec1,data_vec2,data_vec2);
q=[q1(:),q2(:)];
dq=[dq1(:),dq2(:)];
ddq=[ddq1(:),ddq2(:)];

tmp=length(q1(:));
y=zeros(tmp,n);
rng(1);
for i=1:tmp
    
y(i,:)=H(q(i,:)')*ddq(i,:)'+C(dq(i,:)',q)*dq(i,:)'+G(q(i,:)')+F(ddq(i,:)',dq(i,:)',q(i,:)')-(hatH(q(i,:)')*ddq(i,:)'+hatC(dq(i,:)',q)*dq(i,:)'+hatG(q(i,:)'))+randn(n,1)/10;
end
disp('...done');
%% GP1 training
disp('Starting GP1 training');

likfunc = @likGauss;
covfunc = @covSEard;

    
%Input data
training_input=[ddq,dq,q];
training_output1=y(:,1);
training_output2=y(:,2);


%GP Training
clear hyp hyp2 hyp3 hyp4 hyp5 hyp6
hyp.cov = ones(7,1); hyp.lik = log(1);
hyp = minimize(hyp, @gp, -20, @infExact, [], covfunc, likfunc, training_input, training_output1);
hyp.sd=sqrt(exp(2*hyp.cov(end)));
hyp.ell=exp(hyp.cov(1:end-1));
hyp.K=inv(kernel_se(training_input',training_input', hyp.sd,hyp.ell)+eye(length(training_input))*exp(2*hyp.lik));
hyp.Ky=hyp.K*training_output1;

hyp2.cov = ones(7,1); hyp2.lik = log(1);
hyp2 = minimize(hyp2, @gp, -20, @infExact, [], covfunc, likfunc, training_input, training_output2);
hyp2.sd=sqrt(exp(2*hyp2.cov(end)));
hyp2.ell=exp(hyp2.cov(1:end-1));
hyp2.K=inv(kernel_se(training_input',training_input', hyp2.sd,hyp2.ell)+eye(length(training_input))*exp(2*hyp2.lik));
hyp2.Ky=hyp2.K*training_output2;

f=@(mod,t,q,dq,ddq) [predict_mod(mod,hyp,training_input',[ddq;dq;q]),predict_mod(mod,hyp2,training_input',[ddq;dq;q])];

%GP Test
[tmp,~] = gp(hyp, @infExact, [], covfunc, likfunc, training_input, training_output1, training_input);

disp('...done');

figure(2);
plot(abs(tmp-training_output1));
title('Training error');
xlabel('time');
ylabel('GPR mean - Y');
drawnow();


disp('Starting GP2 training');

hyp3.lik=hyp2.lik;
hyp3.sd=sqrt(exp(2*hyp2.cov(end)));
hyp3.ell=hyp2.ell(n+1:end);
hyp3.K=inv(kernel_se(training_input(:,n+1:end)',training_input(:,n+1:end)', hyp3.sd,hyp3.ell)+eye(length(training_input(:,n+1:end)))*exp(2*hyp3.lik));
hyp3.Ky=hyp3.K*training_output1;

hyp4.lik=hyp2.lik;
hyp4.sd=sqrt(exp(2*hyp2.cov(end)));
hyp4.ell=hyp2.ell(n+1:end);
hyp4.K=inv(kernel_se(training_input(:,n+1:end)',training_input(:,n+1:end)', hyp4.sd,hyp4.ell)+eye(length(training_input(:,n+1:end)))*exp(2*hyp4.lik));
hyp4.Ky=hyp4.K*training_output2;

hyp5.lik=hyp2.lik;
hyp5.sd=sqrt(exp(2*hyp2.cov(end)));
hyp5.ell=hyp2.ell(2*n+1);
hyp5.K=inv(kernel_se(training_input(:,2*n+1)',training_input(:,2*n+1)', hyp5.sd,hyp5.ell)+eye(length(training_input(:,2*n+1)))*exp(2*hyp5.lik));
hyp5.Ky=hyp5.K*training_output1;

hyp6.lik=hyp2.lik;
hyp6.sd=sqrt(exp(2*hyp2.cov(end)));
hyp6.ell=hyp2.ell(end);
hyp6.K=inv(kernel_se(training_input(:,end)',training_input(:,end)', hyp6.sd,hyp6.ell)+eye(length(training_input(:,end)))*exp(2*hyp6.lik));
hyp6.Ky=hyp6.K*training_output2;

disp('...done');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control with GP
Kpf=@(q) Kp+400*([predict_mod(3,hyp5,training_input(:,2*n+1)',q(1)),0;0,predict_mod(3,hyp6,training_input(:,end)',q(2))])-diag([3,3]);
Kdf=@(dq,q) Kd+400*([predict_mod(3,hyp3,training_input(:,n+1:end)',[dq;q]),0;0,predict_mod(3,hyp4,training_input(:,n+1:end)',[dq;q])])-diag([4,4]);
u_gp_mean=@(t,ddq,dq,q) hatH(q)*ddqd(t)+hatC(dq,q)*dqd(t)+hatG(q)+f(2,t,q,dq,ddq)'-Kpf(q)*(q-qd(t))-Kdf(dq,q)*(dq-dqd(t));

%% Simulation ODE

disp('Starting simulation');

[t,x] = ode15i(@(t,x,xp) sde_robot(t,x,xp,n,u_gp_mean,H,C,G,F), tspan,x0,xp0,opts);

disp('...done');
%% Plot

figure(3);
clf
title('Control with CTC-GPR')
subplot(2,1,1)
plot(t,qd_vec(:,1),'--',t,dqd_vec(:,1),'--');
hold on
plot(t,x(:,[1 3]));
legend('dq1','dqd1','q1','qd1');
xlabel('time');
ylabel('States');

subplot(2,1,2)
plot(t,qd_vec(:,2),'--',t,dqd_vec(:,2),'--');
hold on
plot(t,x(:,[2 4]));
legend('dq2','dqd2','q2','qd2');
xlabel('time');
ylabel('States');

drawnow();



% xy plot
normgain=zeros(length(x),1);
normgainKd=zeros(length(x),1);
normgainKp=zeros(length(x),1);
normgainKd1=zeros(length(x),1);
normgainKp1=zeros(length(x),1);
normgainKd2=zeros(length(x),1);
normgainKp2=zeros(length(x),1);
for i=1:length(x)
    temp=t(i);
    q=x(i,1:2)';
    dq=x(i,3:4)';
    tmp=Kdf(dq,q);
    tmp1=Kpf(q);
    
    normgainKd1(i)=tmp(1,1);
    normgainKp1(i)=tmp1(1,1);
    normgainKd2(i)=tmp(2,2);
    normgainKp2(i)=tmp1(2,2);
    normgain(i)=normgainKp1(i)+normgainKd1(i);
end

figure(4);
clf
plot(t,normgainKp1,t,normgainKd1,t,Kp(1,1)*ones(length(t),1),'--',t,Kd(1,1)*ones(length(t),1),'--');
title('Adapted feedback gains')
xlabel('time');
ylabel('gains');
legend({'CTC-GPR Kp1','CTC-GPR Kd1','CTC Kp','CTC Kd'});
%
figure(5);
clf
plot(training_input(:,5),training_input(:,3),'+','color',[0 0.7 0.3],'markersize',15,'linewidth',2)
hold on
plot(qd_vec(1:130,1),dqd_vec(1:130,1),'--','linewidth',1);
plot3(x_ct(:,1),x_ct(:,n+1),ones(length(x_ct(:,1)),1)*norm(Kp)+norm(Kd));
c=normgain;
colormap('jet');
cmap = colormap;
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
plot3(x(:,1),x(:,n+1),normgain,'linestyle','none','HandleVisibility','off')
hold on;
for k = 1:(length(x)-1)
    line(x(k:k+1,1),x(k:k+1,n+1),normgain(k:k+1),'color',cmap(c(k),:),'linewidth',2)
end
caxis([ min(c) , max(c)])

legend('Training data','desired trajectory','classical CTC','CTC with GPR');
title('Closed loop error for first joint')
xlabel('$q_1$','interpreter','latex');
ylabel('$\dot{q}_1$','interpreter','latex');
drawnow();
%% Control with GP constant gains
u_gp_mean_const=@(t,ddq,dq,q) hatH(q)*ddqd(t)+hatC(dq,q)*dqd(t)+hatG(q)+f(2,t,q,dq,ddq)'-[min(normgainKp1),0;0,min(normgainKp2)]*(q-qd(t))-[min(normgainKd1),0;0,min(normgainKd2)]*(dq-dqd(t));

%% Simulation ODE

disp('Starting simulation');

[t,x_ctc] = ode15i(@(t,x,xp) sde_robot(t,x,xp,n,u_gp_mean_const,H,C,G,F), tspan,x0,xp0,opts);
    
disp('...done');

%% tracking error
figure(10)
clf
plot(t,abs(qd_vec(:,1)-x_ct(:,1)),t,abs(qd_vec(:,1)-x_ctc(:,1)),t,abs(qd_vec(:,1)-x(:,1)));
legend('CTC','CTC-GPR const gains','CTC-GPR adapted gains');
title('Tracking error');
xlabel('time');
ylabel('States');