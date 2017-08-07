// genesis
// cell parameter file for the 1991 Traub CA1 hippocampal cell
// "phi" parameter reduced by e-3
*cartesian
*relative

*set_global RM 1.0	//ohm*m^2
*set_global RA 1.0	//ohm*m
*set_global CM 0.03     //F/m^2
*set_global EREST_ACT	-0.06	// volts

// The format for each compartment parameter line is :
// name  parent  x       y       z       d       ch      dens ...
// For channels, "dens" =  maximum conductance per unit area of compartment

// First compartment is soma and then it is divided into branches
soma	none    20      0       0       10
soma_2  soma    10      0       0       10
soma_3  soma_2  10      0       0       10
soma_4  soma_3  10      0       0       10
soma_5  soma_4  10      0       0       10
soma_6  soma_5  10      0       0       10
soma_7  soma_6  10      0       0       10
soma_8  soma_7  20      0       0       10
