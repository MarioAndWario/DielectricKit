/*------------------------------------------------------------------------------

Header for routines common to Visual/surface.cpp and MeanField/ICM/icm.cpp
  which are in wfn_utils.cpp

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
extern "C" {
#include "periodic_table.h"
}

#define MAX(aa,bb) ((aa) > (bb) ? (aa) : (bb))
#define MAXCHAR 256

typedef struct {
   double xx;
   double yy;
   double zz;
} CARTESIAN;

extern const double EPS9;
extern const double INF9;
extern const double BOHR;

/* AAH! global variables! --DAS */
extern int *as;
extern CARTESIAN *ap;

extern const int vertex[8][3];
extern const ELEMENT periodic_table[119];

struct scalar_field
{
   int ni, nj, nk;    // grid size
   CARTESIAN origin;
   CARTESIAN stepv[3];
   double ***scalar;
};

void lowercase(char *ss);
void erasechar(char *ss, char cc);
void parse(char *s1, char *s2, char *s3);
double*** allocate_scalar(int ni, int nj, int nk);
void deallocate_scalar(double ***myscalar, int ni, int nj);
int cub_read(char *ifn, int *na, CARTESIAN *sfv, scalar_field & sf);
int xsf_read(char *ifn, int *na, CARTESIAN *sfv, scalar_field & sf);
void cell_set(CARTESIAN sfo, CARTESIAN *sfv, bool uc, CARTESIAN *uco, CARTESIAN *ucv, int ucf, bool sc, CARTESIAN *sco, CARTESIAN *scv, int scf);
void normal_make(CARTESIAN *vector, CARTESIAN *normal, double *weight);
bool box_check(int ii, int jj, int kk, CARTESIAN *gso, CARTESIAN *sfs, CARTESIAN *sco, CARTESIAN *normal, int *nn);
int scalar_clone(CARTESIAN uco, CARTESIAN *ucv, CARTESIAN sco, CARTESIAN *scv, scalar_field & sf);
double inversion(int nn, double *xx, double *yy, double zz);
double isovalue_scale(int power, double isovalue, const scalar_field & sf);
int double_to_int(double dd);

