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

disp('Compiling mex functions');

cd two_link
mex kernel_se.c
cd ..

cd elsystem
mex kernel_se.c
cd ..

disp('The following functions require the GPML toolbox (http://www.gaussianprocess.org/gpml/code/matlab/doc/)');
disp('Please start two_link\start.m for an example with a robot manipulator');
disp('Please start elsystem\start.m for an example with randomly generated systems');