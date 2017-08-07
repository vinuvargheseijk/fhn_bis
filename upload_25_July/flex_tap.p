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
soma_3  soma_2  1      0       0       1.07
soma_4  soma_3  1      0       0       1.04
soma_5  soma_4  1      0       0       1.01
soma_6  soma_5  1      0       0       0.98
soma_7  soma_6  1      0       0       0.95
soma_8  soma_7  1      0       0       0.92
soma_9  soma_8  1      0       0       0.89
soma_10  soma_9 1      0       0       0.86
soma_11  soma_10 1     0       0       0.83
soma_12  soma_11 1     0       0       0.80
soma_13  soma_12 1     0       0       0.77
soma_14  soma_13 1     0       0       0.74
soma_15  soma_14  1      0      0       0.71
soma_16  soma_15 1     0       0       0.68
soma_17  soma_16  1      0      0       0.65

spine_soma_2_neck       soma_2                  0  0  1.4  0.28
spine_soma_2_head       spine_soma_2_neck       0  0  0.7  0.7

spine_soma_3_neck       soma_3                  0  0  1.4  0.28
spine_soma_3_head       spine_soma_3_neck       0  0  0.7  0.7

spine_soma_4_neck       soma_4                  0  0   1.4  0.28
spine_soma_4_head       spine_soma_4_neck       0  0   0.7  0.7

spine_soma_5_neck       soma_5                  0  0   1.4  0.28
spine_soma_5_head       spine_soma_5_neck       0  0   0.7  0.7

spine_soma_6_neck       soma_6                  0  0   1.4  0.28
spine_soma_6_head       spine_soma_6_neck       0  0   0.7  0.7

spine_soma_7_neck       soma_7                  0  0  1.4  0.28
spine_soma_7_head       spine_soma_7_neck       0  0  0.7  0.7

spine_soma_8_neck       soma_8                  0  0  1.4  0.28
spine_soma_8_head       spine_soma_8_neck       0  0  0.7  0.7

spine_soma_9_neck       soma_9                  0  0   1.4  0.28
spine_soma_9_head       spine_soma_9_neck       0  0   0.7  0.7

spine_soma_10_neck       soma_10                 0  0   1.4  0.28
spine_soma_10_head       spine_soma_10_neck       0  0   0.7  0.7

spine_soma_11_neck       soma_11                 0  0   1.4  0.28
spine_soma_11_head       spine_soma_11_neck      0  0   0.7  0.7

spine_soma_12_neck       soma_12                 0  0   1.4  0.28
spine_soma_12_head       spine_soma_12_neck      0  0   0.7  0.7

spine_soma_13_neck       soma_13                 0  0   1.4  0.28
spine_soma_13_head       spine_soma_13_neck      0  0   0.7  0.7

spine_soma_14_neck       soma_14                 0  0   1.4  0.28
spine_soma_14_head       spine_soma_14_neck      0  0   0.7  0.7

spine_soma_15_neck       soma_15                 0  0   1.4  0.28
spine_soma_15_head       spine_soma_15_neck      0  0   0.7  0.7

spine_soma_16_neck       soma_16                 0  0   1.4  0.28
spine_soma_16_head       spine_soma_16_neck      0  0   0.7  0.7
