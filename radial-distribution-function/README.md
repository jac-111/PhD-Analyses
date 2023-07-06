# Radial Distribution Function (G(r))

The radial distribution function is a measure to understand the packing and order within a molecular system. It provides the probability for an atom to appear at a given distance of separation from the subject atom. This is achieved by first selecting an atom and constructing a sphere of radius ’r’. A second sphere is then constructed of radius r+∆r. All atoms of interest within the interstitial region between the two spheres are then counted. This is iterated over a range of values to sample the surrounding area. Values are normalised by shell volume and averaged over every configuration that forms part of the trajectory file. G(r) can be represented with the following equation:

$$G(r) = \frac{1}{\rho 4 \pi r^{2} \delta r} \sum_{r_{0}=1}^{N-1} \sum_{i=2}^{N} \delta (r=|r_{0}-r_{i}|)$$

where $\rho$ is the density for the system at a shell radius of r with a shell thickness of dr. The below image shows a 2D array of atoms with the sampling section of the G(r) analysis being carried out on the orange atom with an idealised plotted G(r) given that the atoms are constrained (right). 

![Gr_graph](https://github.com/jac-111/PhD-Analyses/blob/main/images/Gr_graph.png?raw=true)