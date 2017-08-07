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
soma	none    1       0       0      1.1
soma_2  soma    1      0       0       1.1
soma_3  soma_2  1      0       0       1.09
soma_4  soma_3  1      0       0       1.08
soma_5  soma_4  1      0       0       1.07
soma_6  soma_5  1      0       0       1.06
soma_7  soma_6  1      0       0       1.05
soma_8  soma_7  1      0       0       1.04
soma_9  soma_8  1      0       0       1.03
soma_10  soma_9 1      0       0       1.02
soma_11  soma_10  1    0       0       1.01
soma_12  soma_11  1    0       0       1.00
soma_13  soma_12  1    0       0       0.99
soma_14  soma_13  1    0       0       0.98
soma_15  soma_14  1    0       0       0.97
soma_16  soma_15  1    0       0       0.96
soma_17  soma_16  1    0       0       0.95
soma_18  soma_17  1    0       0       0.94
soma_19  soma_18  1    0       0       0.93
soma_20  soma_19  1    0       0       0.92
soma_21  soma_20  1    0       0       0.91
soma_22  soma_21  1    0       0       0.91
