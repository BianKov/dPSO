# *d*PSO model
Code for generating hyperbolic networks of **homogeneous angular node distribution** with the *d*-dimensional popularity-similarity optimisation model.  
The functions calculating the cutoff distance in the connection probability function are written in Matlab, the rest of the code is written in Python 3.


# Three-dimensional nPSO model
Code for generating networks based on popularity-similarity optimisation in the 3-dimensional hyperbolic space, using a **mixture of von Mises–Fisher distributions as the angular node distribution**, where denser angular regions can serve as built-in communities.  
The functions calculating the cutoff distance in the connection probability function are written in Matlab, the rest of the code is written in Python 3.


# *f*PSO model
Code for generating **2-dimensional** hyperbolic networks treating the multiplying factor *f* of the initial radial coordinates as a model parameter. In contrast to the original 2-dimensional PSO model, in this version degree decay exponents below 2 are also achievable.  
The code is written in Python 3.


# Output files
- The model parameters are saved in two files: one contains only the values while the other also the name of the parameter for each value
- At temperatures 0<*T*, the programs save for each node the cutoff distance in the connection probability. (Note that at *T*=0, the cutoff distance is irrelevant, and therefore, it is not calculated by the program.)
- The edge lists are saved in text files consisting of two or three columns separated by tabs, where the first 2 columns contain the indexes of the connected nodes and the optional third column contains the link weights calculated as 1/(1+hyperbolic distance between the two nodes in question).
- The node coordinates are saved in text files in which the first column contains the node indexes, the second column contains the radial coordinates in the native representation of the hyperbolic space and the remaining columns contain the (*d* number of) Cartesian coordinates that describe the node positions.


# Reference
arXiv:2108.03328 [physics.soc-ph]

For any problem, please contact Bianka Kovács: <bianka.kovacs@ttk.elte.hu>
