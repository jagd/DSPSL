#ifndef __MOM_MESH_H__
#define __MOM_MESH_H__

/*

           The Port
'''''''''''''''''''''''''''''''''''''''''
'                                       '
'     DIELECTRIC1   STRIP0              '
'    ---------------=======---------    '
'    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\    '
'    \\\\\\\\   SUBSTRATE   \\\\\\\\    '
'    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\    '
'    --------========---------------    '
'            STRIP1    DIELECTRIC2      '
'                                       '
'''''''''''''''''''''''''''''''''''''''''

*/



/*  the interval would be: [START, END) */

#define ID_STRIP_START		0
#define ID_STRIP0_START		0
#define ID_STRIP0_END		1
#define ID_STRIP1_START		1
#define ID_STRIP1_END		2
#define ID_STRIP_END		2
#define ID_DIELECTRIC_START	2
#define ID_DIELECTRIC1_START	3
#define ID_DIELECTRIC1_END	4
#define ID_DIELECTRIC2_START	4
#define ID_DIELECTRIC2_END	5
#define ID_DIELECTRIC_END	5
#define ID_MESH_CELLS		5
#define INDEX_SIZE		6


struct Cell_1D {
	double centre;
	double hw; /* half of the width */
};

struct MeshConfig {
	int index[INDEX_SIZE];
	struct Cell_1D *mesh;
};

#endif
