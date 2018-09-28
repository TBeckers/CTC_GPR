# CTC_GPR
Stable Gaussian Process based Tracking Control of Euler-Lagrange Systems

Perfect tracking control for real-world Euler-Lagrange systems is challenging due to uncertainties in the system model and external disturbances. The magnitude of the tracking error can be reduced either by increasing the feedback gains or improving the model of the system. The latter is clearly preferable as it allows to maintain good tracking performance at low feedback gains. However, accurate models are often difficult to obtain.

In this article, we address the problem of high-performance tracking control for unknown Euler-Lagrange systems. In particular, we employ Gaussian Process Regression (GPR) to obtain a data-driven 		model that is used for the feed-forward compensation of unknown dynamics of the system. Beneficially, GPR provides not only an estimate of the uncertainties, but naturally provides a measure of model confidence depending on the distance to training points. Accordingly, the feedback gain is adapted based on the model fidelity allowing low feedback gains in state space regions of high model confidence. Additionally, we study the stability of GP-based tracking control for Euler-Lagrange systems. The proposed confidence-adaptive feedback control law guarantees a globally bounded tracking error with a specific probability. Simulation studies illustrate the results and demonstrate the superiority over state of the art tracking control approaches.

Please start the m-file "start.m"

Tested with Matlab R2018a / Microsoft Windows 10
