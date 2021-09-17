/******************************************************************************
*
*  spglib_f_meta           Originally by FHJ     Last Modified 11/08/2012 (FHJ)
*
*
*  SHORT
*    BerkeleyGW meta wrapper for spglib.
*
*  LONG
*    spglib`s Fortran wrappers assume the Fortran compiler will append an
*    underscore after each symbol`s name, which is not true for all compilers,
*    such as XL. This can lead to terrible results, specially since the symbols
*    without underscores are also available in the library, but take different
*    arguments. This file fixes this problem by exporting symbols with both
*    "_f" and "_f_" appended to the names, so that they should work regardless
*    of the name mangling.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/


/* Original symbol from spglib_f.c */
extern void spg_get_multiplicity_( int *size,
			    double lattice[3][3],
			    double position[][3],
			    int types[],
			    int *num_atom,
			    double *symprec );
/* Wrapped symbol without underscore */
void spg_get_multiplicity_f( int *size,
			    double lattice[3][3],
			    double position[][3],
			    int types[],
			    int *num_atom,
			    double *symprec ) {
    spg_get_multiplicity_(size,
			   lattice,
			   position,
			   types,
			   num_atom,
			   symprec );
}
/* Wrapped symbol with underscore */
void spg_get_multiplicity_f_( int *size,
			    double lattice[3][3],
			    double position[][3],
			    int types[],
			    int *num_atom,
			    double *symprec ) {
    spg_get_multiplicity_(size,
			   lattice,
			   position,
			   types,
			   num_atom,
			   symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_get_max_multiplicity_( int *size,
				double lattice[3][3],
				double position[][3],
				int types[],
				int *num_atom,
				double *symprec );
/* Wrapped symbol without underscore */
void spg_get_max_multiplicity_f( int *size,
				double lattice[3][3],
				double position[][3],
				int types[],
				int *num_atom,
				double *symprec ) {
    spg_get_max_multiplicity_(size,
			lattice,
			position,
			types,
			num_atom,
			symprec );
}
/* Wrapped symbol with underscore */
void spg_get_max_multiplicity_f_( int *size,
				double lattice[3][3],
				double position[][3],
				int types[],
				int *num_atom,
				double *symprec ) {
    spg_get_max_multiplicity_(size,
			lattice,
			position,
			types,
			num_atom,
			symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_get_symmetry_( int *nsym,
			int rot[][3][3],
			double trans[][3],
			int *size,
			double lattice[3][3],
			double position[][3],
			int types[],
			int *num_atom,
			double *symprec );
/* Wrapped symbol without underscore */
void spg_get_symmetry_f( int *nsym,
			int rot[][3][3],
			double trans[][3],
			int *size,
			double lattice[3][3],
			double position[][3],
			int types[],
			int *num_atom,
			double *symprec ) {
    spg_get_symmetry_(nsym,
		rot,
		trans,
		size,
		lattice,
		position,
		types,
		num_atom,
		symprec );
}
/* Wrapped symbol with underscore */
void spg_get_symmetry_f_( int *nsym,
			int rot[][3][3],
			double trans[][3],
			int *size,
			double lattice[3][3],
			double position[][3],
			int types[],
			int *num_atom,
			double *symprec ) {
    spg_get_symmetry_(nsym,
		rot,
		trans,
		size,
		lattice,
		position,
		types,
		num_atom,
		symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_get_smallest_lattice_( double smallest_lattice[3][3],
				double lattice[3][3],
				double *symprec );
/* Wrapped symbol without underscore */
void spg_get_smallest_lattice_f( double smallest_lattice[3][3],
				double lattice[3][3],
				double *symprec ) {
    spg_get_smallest_lattice_(smallest_lattice,
			lattice,
			symprec );
}
/* Wrapped symbol with underscore */
void spg_get_smallest_lattice_f_( double smallest_lattice[3][3],
				double lattice[3][3],
				double *symprec ) {
    spg_get_smallest_lattice_(smallest_lattice,
			lattice,
			symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_get_international_( int *spacegroup,
			     char symbol[11],
			     double lattice[3][3],
			     double position[][3],
			     int types[],
			     int *num_atom,
			     double *symprec );
/* Wrapped symbol without underscore */
void spg_get_international_f( int *spacegroup,
			     char symbol[11],
			     double lattice[3][3],
			     double position[][3],
			     int types[],
			     int *num_atom,
			     double *symprec ) {
    spg_get_international_(spacegroup,
			    symbol,
			    lattice,
			    position,
			    types,
			    num_atom,
			    symprec );
}
/* Wrapped symbol with underscore */
void spg_get_international_f_( int *spacegroup,
			     char symbol[11],
			     double lattice[3][3],
			     double position[][3],
			     int types[],
			     int *num_atom,
			     double *symprec ) {
    spg_get_international_(spacegroup,
			    symbol,
			    lattice,
			    position,
			    types,
			    num_atom,
			    symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_refine_cell_( double lattice[3][3],
		       double position[][3],
		       int types[],
		       int *num_atom,
		       double *symprec );
/* Wrapped symbol without underscore */
void spg_refine_cell_f( double lattice[3][3],
		       double position[][3],
		       int types[],
		       int *num_atom,
		       double *symprec ) {
    spg_refine_cell_(lattice,
		      position,
		      types,
		      num_atom,
		      symprec );
}
/* Wrapped symbol with underscore */
void spg_refine_cell_f_( double lattice[3][3],
		       double position[][3],
		       int types[],
		       int *num_atom,
		       double *symprec ) {
    spg_refine_cell_(lattice,
		      position,
		      types,
		      num_atom,
		      symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_get_schoenflies_( int *spacegroup,
			   char symbol[10],
			   double lattice[3][3],
			   double position[][3],
			   int types[],
			   int *num_atom,
			   double *symprec );
/* Wrapped symbol without underscore */
void spg_get_schoenflies_f( int *spacegroup,
			   char symbol[10],
			   double lattice[3][3],
			   double position[][3],
			   int types[],
			   int *num_atom,
			   double *symprec ) {
    spg_get_schoenflies_(spacegroup,
			  symbol,
			  lattice,
			  position,
			  types,
			  num_atom,
			  symprec );
}
/* Wrapped symbol with underscore */
void spg_get_schoenflies_f_( int *spacegroup,
			   char symbol[10],
			   double lattice[3][3],
			   double position[][3],
			   int types[],
			   int *num_atom,
			   double *symprec ) {
    spg_get_schoenflies_(spacegroup,
			  symbol,
			  lattice,
			  position,
			  types,
			  num_atom,
			  symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_find_primitive_( double lattice[3][3],
			  double position[][3],
			  int types[],
			  int *num_atom,
			  double *symprec );
/* Wrapped symbol without underscore */
void spg_find_primitive_f( double lattice[3][3],
			  double position[][3],
			  int types[],
			  int *num_atom,
			  double *symprec ) {
    spg_find_primitive_(lattice,
			 position,
			 types,
			 num_atom,
			 symprec );
}
/* Wrapped symbol with underscore */
void spg_find_primitive_f_( double lattice[3][3],
			  double position[][3],
			  int types[],
			  int *num_atom,
			  double *symprec ) {
    spg_find_primitive_(lattice,
			 position,
			 types,
			 num_atom,
			 symprec );
}

/* Original symbol from spglib_f.c */
extern void spg_get_ir_reciprocal_mesh_( int *num_ir_grid,
				  int grid_point[][3],
				  int map[],
				  int mesh[3],
				  int is_shift[3],
				  int *is_time_reversal,
				  double lattice[3][3],
				  double position[][3],
				  int types[],
				  int *num_atom,
				  double *symprec );
/* Wrapped symbol without underscore */
void spg_get_ir_reciprocal_mesh_f( int *num_ir_grid,
				  int grid_point[][3],
				  int map[],
				  int mesh[3],
				  int is_shift[3],
				  int *is_time_reversal,
				  double lattice[3][3],
				  double position[][3],
				  int types[],
				  int *num_atom,
				  double *symprec ) {
    spg_get_ir_reciprocal_mesh_(num_ir_grid,
				 grid_point,
				 map,
				 mesh,
				 is_shift,
				 is_time_reversal,
				 lattice,
				 position,
				 types,
				 num_atom,
				 symprec );
}
/* Wrapped symbol with underscore */
void spg_get_ir_reciprocal_mesh_f_( int *num_ir_grid,
				  int grid_point[][3],
				  int map[],
				  int mesh[3],
				  int is_shift[3],
				  int *is_time_reversal,
				  double lattice[3][3],
				  double position[][3],
				  int types[],
				  int *num_atom,
				  double *symprec ) {
    spg_get_ir_reciprocal_mesh_(num_ir_grid,
				 grid_point,
				 map,
				 mesh,
				 is_shift,
				 is_time_reversal,
				 lattice,
				 position,
				 types,
				 num_atom,
				 symprec );
}
