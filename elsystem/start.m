% Copyright (C) Technical University of Munich
% Written by Thomas Beckers, ITR, Technical University of Munich, 2018

clear RMSE_ct RMSE_gpc SNR_ct SNR_gpc MAX_ct MAX_gpc hyp Kp_gpc Kd_gpc

n_trials=5;

tspan=0:0.01:2*pi;

rng(6)
rand_q0_array=2*pi*rand(n_trials,1);

for trial=1:n_trials
disp(['Starting ' num2str(trial)]);

m=1;
k=1;
b=1;

H=@(q) m;
C=@(dq,q) b;
G=@(q) k;

% Desired trajectory
qd=@(t) sin(t);
dqd=@(t) cos(t);
ddqd=@(t) -sin(t);

% Control law
Kp= 100;
Kd= 100;

rand_q0=rand_q0_array(trial);

F=@(ddq,dq,q) 1*dp(rand_q0,dq,q);


u=@(t,ddq,dq,q) H(q)*ddqd(t)+C(dq,q)*dqd(t)+G(q)-Kp*(q-qd(t))-Kd*(dq-dqd(t));

%% Simulation with classical CTC


q0=0;
dq0=1;
ddq0=0;

x0=[dq0;q0;0];
xp0=[ddq0;dq0;0];

[t,x_ct] = ode45(@(t,x) sde_robot(t,x,u,H,C,G,F), tspan,x0);

u_ct= odepassout(@(t,x) sde_robot(t,x,u,H,C,G,F),t,x_ct,1);
% Plot
figure(1);
clf
subplot(2,1,1);
plot(t,dqd(t)','--',t,qd(t)','--');
hold on
plot(t,x_ct);
legend('dqd','qd','dq','q');
title('Control with classical CTC');
xlabel('time');
ylabel('States');

subplot(2,1,2);
plot(t,u_ct);
legend('u');
xlabel('time');
ylabel('control action');
drawnow();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data generation for learning (direct q,dq,ddq)
data_vec=-1:0.1:1;
[ddqt,dqt,qt]=ndgrid(0,data_vec,data_vec);

rng(1);
tmp=length(qt(:));
y=zeros(tmp,1);
for i=1:tmp 
    ddqtrue=ddqt(i);
    dqtrue=dqt(i);
    qtrue=qt(i);
    ddqt(i)=ddqt(i)+randn/100;
    dqt(i)=dqt(i)+randn/100;
    qt(i)=qt(i)+randn/100;
    y(i)=H(qtrue)*ddqtrue+C(dqtrue,qtrue)*dqtrue+G(qtrue)+F(ddqtrue,dqtrue,qtrue)-(H(qt(i))*ddqt(i)+C(dqt(i),qt(i))*dqt(i)+G(qt(i)));
end
y=y+randn(tmp,1)/50;
%% GP training

likfunc = @likGauss;
covfunc = @covSEard;

    
%Input data
training_input=[ddqt(:),dqt(:),qt(:)];

%GP Training
hyp=struct('cov',ones(4,1),'lik',log(1),'sd',[],'ell',[],'K',[],'Ky',[]);

hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, training_input, y);
hyp.sd=sqrt(exp(2*hyp.cov(end)));
hyp.ell=exp(hyp.cov(1:end-1));
hyp.K=inv(kernel_se(training_input',training_input', hyp.sd,hyp.ell)+eye(length(training_input))*exp(2*hyp.lik));
hyp.Ky=hyp.K*y;

hypq=[];
hypq.sd=hyp.sd;
hypq.ell=hyp.ell(3);
hypq.K=inv(kernel_se(training_input(:,3)',training_input(:,3)', hypq.sd,hypq.ell)+eye(length(training_input))*exp(2*hyp.lik));

hypdq=[];
hypdq.sd=hyp.sd;
hypdq.ell=hyp.ell(2:3);
hypdq.K=inv(kernel_se(training_input(:,2:3)',training_input(:,2:3)', hypdq.sd,hypdq.ell)+eye(length(training_input))*exp(2*hyp.lik));

%GP Test
[tmp,~] = gp(hyp, @infExact, [], covfunc, likfunc, training_input, y, training_input);
figure(2);
clf
plot(abs(tmp-y));
title('Training error');
xlabel('time');
ylabel('GPR mean - Y');
drawnow();

f=@(mod,t,ddq,dq,q) predict_mod(mod,hyp,training_input',[ddq;dq;q]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control with CTC-GPR
Kplow=@(q) 10+100*abs(predict_mod(3,hypq,training_input(:,3)',q));
Kdlow=@(dq,q) 10+100*abs(predict_mod(3,hypdq,training_input(:,2:3)',[dq;q]));
u_gp_mean_const=@(t,ddq,dq,q) H(q)*ddqd(t)+C(dq,q)*dqd(t)+G(q)+f(2,t,ddq,dq,q)'-Kplow(q)*(q-qd(t))-Kdlow(dq,q)*(dq-dqd(t));


[t,x_gpc] = ode45(@(t,x) sde_robot(t,x,u_gp_mean_const,H,C,G,F), tspan,x0);
    
u_gpc= odepassout(@(t,x) sde_robot(t,x,u_gp_mean_const,H,C,G,F),t,x_gpc,1);

qd_vec=[dqd(t')',qd(t')'];

% Plot

figure(3);
clf
subplot(2,1,1);
plot(t,qd_vec,'--');
hold on
plot(t,x_gpc);
legend('dqd','qd','dq','q');
title('Control with CTC-GPR')
xlabel('time');
ylabel('States');
subplot(2,1,2);
plot(t,u_gpc);
legend('u');
xlabel('time');
ylabel('control action');
drawnow();

% tracking error
figure(5)
clf
plot(t,rownorm(qd_vec-x_ct(:,1:2)),t,rownorm(qd_vec-x_gpc(:,1:2)));
legend('CTC','CTC-GPR');
title('Tracking error for CTC and CTC-GPR');
xlabel('time');
ylabel('Tracking error');

%% Evaluation

%RMSE_ct(trial) = sqrt(mean(rownorm(qd_vec-x_ct(:,1:2)).^2));
RMSE_ct(trial) = max(rownorm(qd_vec-x_ct(:,1:2)));

%RMSE_gpc(trial) = sqrt(mean(rownorm(qd_vec-x_gpc(:,1:2)).^2));
RMSE_gpc(trial) = max(rownorm(qd_vec-x_gpc(:,1:2)));

SNR_ct(trial)=snr(x_ct(:,1))+snr(x_ct(:,2));

SNR_gpc(trial)=snr(x_gpc(:,1))+snr(x_gpc(:,2));

MAX_ct(trial)=max(abs(u_ct));

MAX_gpc(trial)=max(abs(u_gpc));

Kp_gpc(:,trial)=Kplow(x_gpc(:,2)');

Kd_gpc(:,trial)=Kdlow(x_gpc(:,1)',x_gpc(:,2)');

end
%% Evaluation Trials
figure(10)
clf
subplot(5,1,1)
plot(RMSE_ct,'-+');
hold on;
plot(RMSE_gpc,'-+');
ylabel('RMSE');
legend('CTC','CTC-GPR');
subplot(5,1,2)
plot(SNR_ct,'-+');
hold on;
plot(SNR_gpc,'-+');
ylabel('SNR');
legend('CTC','CTC-GPR');
subplot(5,1,3)
plot(MAX_ct,'-+');
hold on;
plot(MAX_gpc,'-+');
ylabel('Max control action');
legend('CTC','CTC-GPR');
subplot(5,1,4)
plot(repmat(1:n_trials,length(tspan),1),Kp_gpc,'-+');
ylabel('Kp');
subplot(5,1,5)
plot(repmat(1:n_trials,length(tspan),1),Kd_gpc,'-+');
ylabel('Kd');
xlabel('Trials');
drawnow();