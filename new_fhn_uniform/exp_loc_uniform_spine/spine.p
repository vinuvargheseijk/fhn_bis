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
soma	none    20       0       0      10
soma_2  soma    1      0       0       10
soma_3  soma_2  1      0       0       10
soma_4  soma_3  1      0       0       10
soma_5  soma_4  1      0       0       10
soma_6  soma_5  1      0       0       10
soma_7  soma_6  1      0       0       10
soma_8  soma_7  1      0       0       10
soma_9  soma_8  1      0       0       10
soma_10  soma_9 1      0       0       10
soma_11  soma_10  1    0       0       10
soma_12  soma_11  1    0       0       10
soma_13  soma_12  1    0       0       10
soma_14  soma_13  1    0       0       10
soma_15  soma_14  1    0       0       10
soma_16  soma_15  1    0       0       10
soma_17  soma_16  1    0       0       10
soma_18  soma_17  1    0       0       10
soma_19  soma_18  1    0       0       10
soma_20  soma_19  1    0       0       10
soma_21  soma_20  1    0       0       10
soma_22  soma_21  1    0       0       10
soma_23  soma_22  1    0       0       10
soma_24  soma_23  1    0       0       10
soma_25  soma_24  1    0       0       10
soma_26  soma_25  1    0       0       10
soma_27  soma_26  1    0       0       10
soma_28  soma_27  1    0       0       10
soma_29  soma_28  1    0       0       10
soma_30  soma_29  1    0       0       10
soma_31  soma_30  1    0       0       10
soma_32  soma_31  50    0       0       10

spine_soma_neck       soma                    0  0  6.4  0.28
spine_soma_head       spine_soma_neck         0  0  0.7  0.7

spine_soma_2_neck       soma_2                  0  0   6.4  0.28
spine_soma_2_head       spine_soma_2_neck       0  0   0.7  0.7

spine_soma_3_neck       soma_3                  0  0   6.4  0.28
spine_soma_3_head       spine_soma_3_neck       0  0   0.7  0.7

spine_soma_4_neck       soma_4                  0  0   6.4  0.28
spine_soma_4_head       spine_soma_4_neck       0  0   0.7  0.7

spine_soma_5_neck       soma_5                  0  0  6.4  0.28
spine_soma_5_head       spine_soma_5_neck       0  0  0.7  0.7

spine_soma_6_neck       soma_6                  0  0  6.4  0.28
spine_soma_6_head       spine_soma_6_neck       0  0  0.7  0.7

spine_soma_7_neck       soma_7                  0  0   6.4  0.28
spine_soma_7_head       spine_soma_7_neck       0  0   0.7  0.7

spine_soma_8_neck       soma_8                 0  0   6.4  0.28
spine_soma_8_head       spine_soma_8_neck       0  0   0.7  0.7

spine_soma_9_neck       soma_9                 0  0   6.4  0.28
spine_soma_9_head       spine_soma_9_neck      0  0   0.7  0.7

spine_soma_10_neck       soma_10                 0  0   6.4  0.28
spine_soma_10_head       spine_soma_10_neck      0  0   0.7  0.7

spine_soma_11_neck       soma_11                 0  0   6.4  0.28
spine_soma_11_head       spine_soma_11_neck      0  0   0.7  0.7

spine_soma_12_neck       soma_12                 0  0   6.4  0.28
spine_soma_12_head       spine_soma_12_neck      0  0   0.7  0.7

spine_soma_13_neck       soma_13                 0  0   6.4  0.28
spine_soma_13_head       spine_soma_13_neck      0  0   0.7  0.7

spine_soma_14_neck       soma_14                 0  0   6.4  0.28
spine_soma_14_head       spine_soma_14_neck      0  0   0.7  0.7

spine_soma_15_neck       soma_15                 0  0   6.4  0.28
spine_soma_15_head       spine_soma_15_neck      0  0   0.7  0.7

spine_soma_16_neck       soma_16                 0  0   6.4  0.28
spine_soma_16_head       spine_soma_16_neck      0  0   0.7  0.7

spine_soma_17_neck       soma_17                 0  0   6.4  0.28
spine_soma_17_head       spine_soma_17_neck      0  0   0.7  0.7

spine_soma_18_neck       soma_18                 0  0   6.4  0.28
spine_soma_18_head       spine_soma_18_neck      0  0   0.7  0.7

spine_soma_19_neck       soma_19                 0  0   6.4  0.28
spine_soma_19_head       spine_soma_19_neck      0  0   0.7  0.7

spine_soma_20_neck       soma_20                 0  0   6.4  0.28
spine_soma_20_head       spine_soma_20_neck      0  0   0.7  0.7

spine_soma_21_neck       soma_21                 0  0   6.4  0.28
spine_soma_21_head       spine_soma_21_neck      0  0   0.7  0.7

spine_soma_22_neck       soma_22                 0  0   6.4  0.28
spine_soma_22_head       spine_soma_22_neck      0  0   0.7  0.7

spine_soma_23_neck       soma_23                 0  0   6.4  0.28
spine_soma_23_head       spine_soma_23_neck      0  0   0.7  0.7

spine_soma_24_neck       soma_24                 0  0   6.4  0.28
spine_soma_24_head       spine_soma_24_neck      0  0   0.7  0.7

spine_soma_25_neck       soma_25                 0  0   6.4  0.28
spine_soma_25_head       spine_soma_25_neck      0  0   0.7  0.7

spine_soma_26_neck       soma_26                 0  0   6.4    0.28
spine_soma_26_head       spine_soma_26_neck      0  0   0.7    0.7
