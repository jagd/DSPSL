#ifndef __MOM_MESH_H__
#define __MOM_MESH_H__

/* #define MOM_MESH_ENABLE_DEBUG */

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
#define ID_DIELECTRIC0_START	2
#define ID_DIELECTRIC0_END	3
#define ID_DIELECTRIC1_START	3
#define ID_DIELECTRIC1_END	4
#define ID_DIELECTRIC_END	4
#define ID_MESH_CELLS		4
#define INDEX_SIZE		5


/* relative permittivity under this value will be treated as 1.0 */
#define EPS_MIN 1.35

struct Cell_1D {
	double centre;
	double hw; /* half of the width */
};

struct MeshConfig {
	double h;
	double eps_r;
	int index[INDEX_SIZE];
	struct Cell_1D *mesh;
};


struct MeshConfig* mesh_new(
		double w0, /* strip0.width */
		double w1, /* strip1.width */
		double offset, /* offset of the centre of both strip */
		double w_port_ext, /* extend width in both side */
		double h, /* only for determination of the mesh length */
		double eps_r, /* epsilon_r */
		double mesh_step /* if this value not positive,
				    a mesh step will be suggested */
		);

void mesh_free(struct MeshConfig *m);

#endif
