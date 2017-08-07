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
soma	none    30 	0  	0  	10
// dend_1	soma	10	0	0	10
// dend_2  dend_1 	10	0	0	10

//dend_1  soma    10  0  	0  	10

dend_1_1	soma	15	25.98	0	6.3
// dend_1_2	dend_1_1	5	8.66	0	10
// dend_1_3	dend_1_2	5	8.66	0	10

dend_2_1	soma	15	-25.98	0	6.3
// dend_2_2	dend_2_1	5	-8.66	0	10
// dend_2_3	dend_2_2	5	-8.66	0	10
// dend_4_2	dend_4_1	-5	8.66	0	10
// dend_4_3	dend_4_2	-5	8.66	0	10	
// dend_4_4	dend_4_3	-5	8.66	0	10
// // soma	dend_8    0  12.5  0  8.46


// apical_10 	dend_8  0  12  0  3.0
// apical_11  apical_10  0   12  0  3
// apical_12	apical_11  0   12  0  3
// apical_13  	apical_12  0   12  0  3
// apical_14  	apical_13  0   12  0  3
// apical_15	apical_14  0   12  0  3
// apical_16	apical_15  0   12  0  3
// apical_17	apical_16  0   12  0  3
// apical_18	apical_17  0   12  0  3
// apical_19 apical_18	   0  12   0  3

// apical_11_1  apical_10  -6   6  0  3
// apical_11_2  	apical_11_1  -6   6  0  3
// apical_11_3  	apical_11_2  0   8  0  3
// apical_11_4  	apical_11_3  0   8  0  3

// apical_13_1	apical_12  4   4  0  3
// apical_13_2	apical_13_1  4   4  0  3

// apical_14_1	apical_13_2  0   1  0  3
// apical_14_2  	.  0   1  0  3
// apical_14_3  	.  0   1  0  3
// apical_14_4  	.  0   1  0  3
// apical_14_5  	.  0   1  0  3
// apical_14_6  	.  0   1  0  3
// apical_14_7  	.  0   1  0  3
// apical_14_8  	.  0   1  0  3
// apical_14_9  	.  0   1  0  3
// apical_14_10  	.  0   1  0  3
// apical_14_11  	.  0   1 0  3
// apical_14_12  	.  0   1  0  3
// apical_14_13  	.  0   1  0  3
// apical_15_1  	.  0   3  0  3
// apical_15_2  	.  0   6  0  3

