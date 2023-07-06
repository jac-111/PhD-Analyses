## Mean squared displacement (MSD)
MSD is a means to study the average displacement of a particle from
a given starting position over a range of time. This measure enables the observation of
molecular mobility within a system. For the molecular system below, a very steep MSD curve would be produced. 

![molecues_lowD](https://github.com/jac-111/PhD-Analyses/blob/main/images/molecues_lowD.gif?raw=true)

The equation for MSD is:
$$MSD(τ) = \frac{1}{N_{conf}N_{mol}}\sum^{conf}_{t=0}\sum^{conf}_{t=1}|x_{i}(t+\tau) - x_{i}(t)|^{2}$$

where N<sub>conf</sub> is the number of configurations and N<sub>mol</sub> is the number of molecules respectively, x_{i} is the position of particle i, t is the time and τ is a time interval for each measurement. Each value of τ is iterated over for all configurations.