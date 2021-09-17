/*------------------------------------------------------------------------------

Routines common to Visual/surface.cpp and MeanField/ICM/icm.cpp
   written by Georgy Samsonidze (October 2008)
   refactored into common file by DAS

------------------------------------------------------------------------------*/

#include "wfn_utils.h"

const double EPS9 = 1.0E-9;
const double INF9 = 1.0E+9;
const double BOHR = 0.52917721092;

/* AAH! global variables! --DAS */
int *as = NULL;
CARTESIAN *ap = NULL;

const int vertex[8][3] = {
   {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
   {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};

void lowercase(char *ss)
{
   int ii, nn;

   nn = (int)strlen(ss);
   for (ii = 0; ii < nn; ii++)
   {
      if (ss[ii] >= 'A' && ss[ii] <= 'Z')
         ss[ii] = ss[ii] + 32;
   }
}

void erasechar(char *ss, char cc)
{
   int ii, jj, nn;

   nn = (int)strlen(ss);
   for (ii = nn - 1; ii >= 0; ii--)
      if (ss[ii] == cc)
      {
         nn--;
         for (jj = ii; jj < nn; jj++)
            ss[jj] = ss[jj + 1];
         ss[nn] = '\0';
      }
}

void parse(char *s1, char *s2, char *s3)
{
   int ii;
   char s4[MAXCHAR], s5[MAXCHAR], s6[MAXCHAR], s7[MAXCHAR];

   ii = (int)strspn(s1, " \t");
   strncpy(s4, &s1[ii], MAXCHAR - ii);
   s4[MAXCHAR - 1 - ii] = '\0';
   ii = (int)strcspn(s4, " \t");
   strncpy(s5, s4, ii);
   s5[ii] = '\0';
   erasechar(s5, '\n');
   erasechar(s5, '\r');
   strncpy(s6, &s4[ii], MAXCHAR - ii);
   s6[MAXCHAR - 1 - ii] = '\0';
   ii = (int)strspn(s6, " \t");
   strncpy(s7, &s6[ii], MAXCHAR - ii);
   s7[MAXCHAR - 1 - ii] = '\0';
   erasechar(s7, '\n');
   erasechar(s7, '\r');
   strncpy(s2, s5, MAXCHAR);
   s2[MAXCHAR - 1] = '\0';
   strncpy(s3, s7, MAXCHAR);
   s3[MAXCHAR - 1] = '\0';
}

double*** allocate_scalar(int ni, int nj, int nk)
{
   int ii, jj;
   double ***myscalar;

   myscalar = new double**[ni];
   for (ii = 0; ii < ni; ii++)
      myscalar[ii] = new double*[nj];
   for (ii = 0; ii < ni; ii++)
   for (jj = 0; jj < nj; jj++)
      myscalar[ii][jj] = new double[nk];

   return myscalar;
}

void deallocate_scalar(double ***myscalar, int ni, int nj)
{
   int ii, jj;

   for (ii = 0; ii < ni; ii++)
      for (jj = 0; jj < nj; jj++)
         delete [] myscalar[ii][jj];
   for (ii = 0; ii < ni; ii++)
      delete [] myscalar[ii];
   delete [] myscalar;

   return;
}

// Read a scalar field into global variable "scalar" from a CUBE file.
// ifn = filename
// na = number of atoms
// sfv = scalar field lattice vectors
int cub_read(char *ifn, int *na, CARTESIAN *sfv, scalar_field & sf)
{
   int ii, jj, kk, ierr;
   double ac;
   char s1[MAXCHAR], s2[MAXCHAR];
   FILE *hh;

   hh = fopen(ifn, "r");
   if (hh == NULL)
      return -1;

   if (!fgets(s1, MAXCHAR, hh)) return -1;
   if (!fgets(s2, MAXCHAR, hh)) return -1;
   ierr = fscanf(hh, "%i%le%le%le\n", &(*na), &sf.origin.xx, &sf.origin.yy, &sf.origin.zz);
   ierr = fscanf(hh, "%i%le%le%le\n", &(sf.ni), &sf.stepv[0].xx, &sf.stepv[0].yy, &sf.stepv[0].zz);
   ierr = fscanf(hh, "%i%le%le%le\n", &(sf.nj), &sf.stepv[1].xx, &sf.stepv[1].yy, &sf.stepv[1].zz);
   ierr = fscanf(hh, "%i%le%le%le\n", &(sf.nk), &sf.stepv[2].xx, &sf.stepv[2].yy, &sf.stepv[2].zz);

   if (*na < 0 || sf.ni < 1 || sf.nj < 1 || sf.nk < 1)
      return -1;

   as = new int[MAX(1,*na)];
   ap = new CARTESIAN[MAX(1,*na)];

   for (ii = 0; ii < *na; ii++)
      ierr = fscanf(hh, "%i%le%le%le%le\n", &as[ii], &ac, &ap[ii].xx, &ap[ii].yy, &ap[ii].zz);

   for (ii = 0; ii < *na; ii++)
      if (as[ii] < 0 || as[ii] > 118)
         return -1;

   sf.scalar = allocate_scalar(sf.ni, sf.nj, sf.nk);

   for (ii = 0; ii < sf.ni; ii++)
   for (jj = 0; jj < sf.nj; jj++)
   for (kk = 0; kk < sf.nk; kk++)
      ierr = fscanf(hh, "%le", &sf.scalar[ii][jj][kk]);

   ierr = fclose(hh);
   if (ierr != 0)
      return -1;

   sf.origin.xx *= BOHR;
   sf.origin.yy *= BOHR;
   sf.origin.zz *= BOHR;
   for (ii = 0; ii < 3; ii++)
   {
      sf.stepv[ii].xx = sf.stepv[ii].xx * BOHR;
      sf.stepv[ii].yy = sf.stepv[ii].yy * BOHR;
      sf.stepv[ii].zz = sf.stepv[ii].zz * BOHR;
   }
   for (ii = 0; ii < *na; ii++)
   {
      ap[ii].xx = ap[ii].xx * BOHR;
      ap[ii].yy = ap[ii].yy * BOHR;
      ap[ii].zz = ap[ii].zz * BOHR;
   }

   sfv[0].xx = sf.stepv[0].xx * double(sf.ni);
   sfv[0].yy = sf.stepv[0].yy * double(sf.ni);
   sfv[0].zz = sf.stepv[0].zz * double(sf.ni);
   sfv[1].xx = sf.stepv[1].xx * double(sf.nj);
   sfv[1].yy = sf.stepv[1].yy * double(sf.nj);
   sfv[1].zz = sf.stepv[1].zz * double(sf.nj);
   sfv[2].xx = sf.stepv[2].xx * double(sf.nk);
   sfv[2].yy = sf.stepv[2].yy * double(sf.nk);
   sfv[2].zz = sf.stepv[2].zz * double(sf.nk);

   return 0;
}

// Read a scalar field into global variable "scalar" from an XSF file.
// ifn = filename
// na = number of atoms
// sfv = scalar field lattice vectors
int xsf_read(char *ifn, int *na, CARTESIAN *sfv, scalar_field & sff)
{
   int ii, jj, kk, pp, ierr;
   int num_elements = get_num_elements();
   double sf;
   char ss[MAXCHAR];
   char dummy[MAXCHAR];
   char symbol[4];
   FILE *hh;

   hh = fopen(ifn, "r");
   if (hh == NULL)
      return -1;

   do
   {
      if (!fgets(ss, MAXCHAR, hh)) return -1;
      lowercase(ss);
      erasechar(ss, ' ');
   }
   while (strncmp(ss, "primcoord", 9) != 0);

   ierr = fscanf(hh, "%i%i\n", &(*na), &pp);

   if (*na < 0 || pp != 1)
      return -1;

   as = new int[MAX(1,*na)];
   ap = new CARTESIAN[MAX(1,*na)];

   for (ii = 0; ii < *na; ii++)
   {
      // read string to buffer first, then read atomic species and
      // positions from the buffer, skip reading optional atomic forces,
      // if atomic species are given by symbols convert them to numbers
      if (!fgets(ss, MAXCHAR, hh)) return -1;
      ierr = sscanf(ss, "%s%le%le%le", dummy, &ap[ii].xx, &ap[ii].yy, &ap[ii].zz);
      as[ii] = 0;
      ierr = sscanf(dummy, "%i", &as[ii]);
      if (as[ii] == 0)
      {
         lowercase(dummy);
         erasechar(dummy, ' ');
         for (jj = 0; jj < num_elements; jj++)
         {
            strncpy(symbol, periodic_table[jj].symbol, 4);
            symbol[3] = '\0';
            lowercase(symbol);
            if (strlen(symbol) == strlen(dummy) && strncmp(symbol, dummy, 4) == 0)
            {
               as[ii] = jj;
               break;
            }
         }
      }
   }

   for (ii = 0; ii < *na; ii++)
      if (as[ii] < 0 || as[ii] > 118)
         return -1;

   do
   {
      if (!fgets(ss, MAXCHAR, hh)) return -1;
      lowercase(ss);
      erasechar(ss, ' ');
   }
   while (strncmp(ss, "begin_datagrid_3d", 17) != 0 && strncmp(ss, "datagrid_3d", 11) != 0);

   ierr = fscanf(hh, "%i%i%i\n", &(sff.ni), &(sff.nj), &(sff.nk));
   ierr = fscanf(hh, "%le%le%le\n", &sff.origin.xx, &sff.origin.yy, &sff.origin.zz);
   ierr = fscanf(hh, "%le%le%le\n", &sfv[0].xx, &sfv[0].yy, &sfv[0].zz);
   ierr = fscanf(hh, "%le%le%le\n", &sfv[1].xx, &sfv[1].yy, &sfv[1].zz);
   ierr = fscanf(hh, "%le%le%le\n", &sfv[2].xx, &sfv[2].yy, &sfv[2].zz);

   // general to periodic grid conversion
   // http://www.xcrysden.org/doc/XSF.html
   sff.ni -= 1;
   sff.nj -= 1;
   sff.nk -= 1;

   if (sff.ni < 1 || sff.nj < 1 || sff.nk < 1)
      return -1;

   sff.scalar = allocate_scalar(sff.ni, sff.nj, sff.nk);

   for (kk = 0; kk <= sff.nk; kk++)
   for (jj = 0; jj <= sff.nj; jj++)
   for (ii = 0; ii <= sff.ni; ii++)
   {
      ierr = fscanf(hh, "%le", &sf);
      if (ii < sff.ni && jj < sff.nj && kk < sff.nk)
         sff.scalar[ii][jj][kk] = sf;
   }

   ierr = fclose(hh);
   if (ierr != 0)
      return -1;

   sff.stepv[0].xx = sfv[0].xx / double(sff.ni);
   sff.stepv[0].yy = sfv[0].yy / double(sff.ni);
   sff.stepv[0].zz = sfv[0].zz / double(sff.ni);
   sff.stepv[1].xx = sfv[1].xx / double(sff.nj);
   sff.stepv[1].yy = sfv[1].yy / double(sff.nj);
   sff.stepv[1].zz = sfv[1].zz / double(sff.nj);
   sff.stepv[2].xx = sfv[2].xx / double(sff.nk);
   sff.stepv[2].yy = sfv[2].yy / double(sff.nk);
   sff.stepv[2].zz = sfv[2].zz / double(sff.nk);

   return 0;
}

// Perform unit conversions and convert to cartesian from lattice coordinates if necessary
// sfo = scalar field origin
// sfv = scalar field lattice vectors
// uc = unit cell. If true, uco/ucv were read from parameter file.
// uco = unit cell origin.
// ucv = unit cell lattice vectors.
// ucf = unit cell format
// sc = supercell. If true, supercell is being constructed.
// sco = supercell origin.
// scv = supercell lattice vectors.
// scf = supercell format.
void cell_set(CARTESIAN sfo, CARTESIAN *sfv, bool uc, CARTESIAN *uco, CARTESIAN *ucv, int ucf, bool sc, CARTESIAN *sco, CARTESIAN *scv, int scf)
{
   int ii;
   double c1, c2, c3;

   if (uc)
   {
      if (ucf == 0)
      {
         uco->xx = uco->xx * BOHR;
         uco->yy = uco->yy * BOHR;
         uco->zz = uco->zz * BOHR;
         for (ii = 0; ii < 3; ii++)
         {
            ucv[ii].xx = ucv[ii].xx * BOHR;
            ucv[ii].yy = ucv[ii].yy * BOHR;
            ucv[ii].zz = ucv[ii].zz * BOHR;
         }
      }
      else if (ucf == 2)
      {
         c1 = uco->xx;
         c2 = uco->yy;
         c3 = uco->zz;
         uco->xx = c1 * sfv[0].xx + c2 * sfv[1].xx + c3 * sfv[2].xx;
         uco->yy = c1 * sfv[0].yy + c2 * sfv[1].yy + c3 * sfv[2].yy;
         uco->zz = c1 * sfv[0].zz + c2 * sfv[1].zz + c3 * sfv[2].zz;
         for (ii = 0; ii < 3; ii++)
         {
            c1 = ucv[ii].xx;
            c2 = ucv[ii].yy;
            c3 = ucv[ii].zz;
            ucv[ii].xx = c1 * sfv[0].xx + c2 * sfv[1].xx + c3 * sfv[2].xx;
            ucv[ii].yy = c1 * sfv[0].yy + c2 * sfv[1].yy + c3 * sfv[2].yy;
            ucv[ii].zz = c1 * sfv[0].zz + c2 * sfv[1].zz + c3 * sfv[2].zz;
         }
      }
   }
   else
   {
      uco->xx = sfo.xx;
      uco->yy = sfo.yy;
      uco->zz = sfo.zz;
      for (ii = 0; ii < 3; ii++)
      {
         ucv[ii].xx = sfv[ii].xx;
         ucv[ii].yy = sfv[ii].yy;
         ucv[ii].zz = sfv[ii].zz;
      }
   }

   if (sc)
   {
      if (scf == 0)
      {
         sco->xx = sco->xx * BOHR;
         sco->yy = sco->yy * BOHR;
         sco->zz = sco->zz * BOHR;
         for (ii = 0; ii < 3; ii++)
         {
            scv[ii].xx = scv[ii].xx * BOHR;
            scv[ii].yy = scv[ii].yy * BOHR;
            scv[ii].zz = scv[ii].zz * BOHR;
         }
      }
      else if (scf == 2)
      {
         c1 = sco->xx;
         c2 = sco->yy;
         c3 = sco->zz;
         sco->xx = c1 * ucv[0].xx + c2 * ucv[1].xx + c3 * ucv[2].xx;
         sco->yy = c1 * ucv[0].yy + c2 * ucv[1].yy + c3 * ucv[2].yy;
         sco->zz = c1 * ucv[0].zz + c2 * ucv[1].zz + c3 * ucv[2].zz;
         for (ii = 0; ii < 3; ii++)
         {
            c1 = scv[ii].xx;
            c2 = scv[ii].yy;
            c3 = scv[ii].zz;
            scv[ii].xx = c1 * ucv[0].xx + c2 * ucv[1].xx + c3 * sfv[2].xx;
            scv[ii].yy = c1 * ucv[0].yy + c2 * ucv[1].yy + c3 * sfv[2].yy;
            scv[ii].zz = c1 * ucv[0].zz + c2 * sfv[1].zz + c3 * sfv[2].zz;
         }
      }
   }
   else
   {
      sco->xx = uco->xx;
      sco->yy = uco->yy;
      sco->zz = uco->zz;
      for (ii = 0; ii < 3; ii++)
      {
         scv[ii].xx = ucv[ii].xx;
         scv[ii].yy = ucv[ii].yy;
         scv[ii].zz = ucv[ii].zz;
      }
   }
}

void normal_make(CARTESIAN *vector, CARTESIAN *normal, double *weight)
{
   int ii;

   normal[0].xx = vector[1].yy * vector[2].zz - vector[1].zz * vector[2].yy;
   normal[0].yy = vector[1].zz * vector[2].xx - vector[1].xx * vector[2].zz;
   normal[0].zz = vector[1].xx * vector[2].yy - vector[1].yy * vector[2].xx;
   normal[1].xx = vector[2].yy * vector[0].zz - vector[2].zz * vector[0].yy;
   normal[1].yy = vector[2].zz * vector[0].xx - vector[2].xx * vector[0].zz;
   normal[1].zz = vector[2].xx * vector[0].yy - vector[2].yy * vector[0].xx;
   normal[2].xx = vector[0].yy * vector[1].zz - vector[0].zz * vector[1].yy;
   normal[2].yy = vector[0].zz * vector[1].xx - vector[0].xx * vector[1].zz;
   normal[2].zz = vector[0].xx * vector[1].yy - vector[0].yy * vector[1].xx;

   for (ii = 0; ii < 3; ii ++)
      weight[ii] = sqrt(normal[ii].xx * normal[ii].xx + normal[ii].yy * normal[ii].yy + normal[ii].zz * normal[ii].zz);
   for (ii = 0; ii < 3; ii ++)
      if (weight[ii] < EPS9)
         weight[ii] = 1.0;
   for (ii = 0; ii < 3; ii ++)
   {
      normal[ii].xx = normal[ii].xx / weight[ii];
      normal[ii].yy = normal[ii].yy / weight[ii];
      normal[ii].zz = normal[ii].zz / weight[ii];
   }

   for (ii = 0; ii < 3; ii ++)
      weight[ii] = vector[ii].xx * normal[ii].xx + vector[ii].yy * normal[ii].yy + vector[ii].zz * normal[ii].zz;
   for (ii = 0; ii < 3; ii ++)
   {
      normal[ii].xx = normal[ii].xx * weight[ii];
      normal[ii].yy = normal[ii].yy * weight[ii];
      normal[ii].zz = normal[ii].zz * weight[ii];
   }

   for (ii = 0; ii < 3; ii ++)
      weight[ii] = normal[ii].xx * normal[ii].xx + normal[ii].yy * normal[ii].yy + normal[ii].zz * normal[ii].zz;
   for (ii = 0; ii < 3; ii ++)
      if (weight[ii] < EPS9)
         weight[ii] = 1.0;
}

bool box_check(int ii, int jj, int kk, CARTESIAN *gso, CARTESIAN *sfs, CARTESIAN *sco, CARTESIAN *normal, int *nn)
{
   bool flag;
   int ll, mm;
   CARTESIAN point;
   double weight[3], proj[3];

   flag = true;

   point.xx = gso[0].xx + sfs[0].xx * double(ii) + sfs[1].xx * double(jj) + sfs[2].xx * double(kk) - sco[0].xx;
   point.yy = gso[0].yy + sfs[0].yy * double(ii) + sfs[1].yy * double(jj) + sfs[2].yy * double(kk) - sco[0].yy;
   point.zz = gso[0].zz + sfs[0].zz * double(ii) + sfs[1].zz * double(jj) + sfs[2].zz * double(kk) - sco[0].zz;

   for (ll = 0; ll < 3; ll++)
      weight[ll] = normal[ll].xx * normal[ll].xx + normal[ll].yy * normal[ll].yy + normal[ll].zz * normal[ll].zz;
   for (ll = 0; ll < 3; ll++)
      proj[ll] = (point.xx * normal[ll].xx + point.yy * normal[ll].yy + point.zz * normal[ll].zz) / weight[ll];
   for (ll = 0; ll < 3; ll++)
   {
// this is similar to double_to_int but we want to convert
// the interval (-1,0] to -1, (0,1] to 0, (1,2] to 1, etc.
      mm = double_to_int(proj[ll] - 0.5);
      if (proj[ll] > 1.0 + EPS9 || proj[ll] < -EPS9)
         flag = false;
      nn[ll] = mm;
   }

   return flag;
}

int scalar_clone(CARTESIAN uco, CARTESIAN *ucv, CARTESIAN sco, CARTESIAN *scv, scalar_field & sf)
{
   int ii, jj, kk, ll, di, dj, dk, ui, uj, uk, sni, snj, snk;
   int nn[3], nmin[3], nmax[3], nuc[3], nsc[3];
   double sfsw[3], ucvw[3], scvw[3], pp;
   CARTESIAN sfsn[3], ucvn[3], scvn[3], ssfo, sccenter, sccorner[8], point;
   double ***scscalar = NULL;

   normal_make(sf.stepv, sfsn, sfsw);
   normal_make(ucv, ucvn, ucvw);
   normal_make(scv, scvn, scvw);

   sccenter.xx = sco.xx + scv[0].xx * 0.5 + scv[1].xx * 0.5 + scv[2].xx * 0.5;
   sccenter.yy = sco.yy + scv[0].yy * 0.5 + scv[1].yy * 0.5 + scv[2].yy * 0.5;
   sccenter.zz = sco.zz + scv[0].zz * 0.5 + scv[1].zz * 0.5 + scv[2].zz * 0.5;

   for (ii = 0; ii < 8; ii++)
   {
      sccorner[ii].xx = sccenter.xx;
      sccorner[ii].yy = sccenter.yy;
      sccorner[ii].zz = sccenter.zz;
      for (jj = 0; jj < 3; jj++)
      {
         sccorner[ii].xx += scv[jj].xx * (double(vertex[ii][jj]) - 0.5);
         sccorner[ii].yy += scv[jj].yy * (double(vertex[ii][jj]) - 0.5);
         sccorner[ii].zz += scv[jj].zz * (double(vertex[ii][jj]) - 0.5);
      }
   }

   for (ii = 0; ii < 3; ii++)
   {
      nmin[ii] = INT_MAX;
      nmax[ii] = INT_MIN;
      for (jj = 0; jj < 8; jj++)
      {
         pp = ((sccorner[jj].xx - sf.origin.xx) * sfsn[ii].xx + (sccorner[jj].yy - sf.origin.yy) * sfsn[ii].yy + (sccorner[jj].zz - sf.origin.zz) * sfsn[ii].zz) / sfsw[ii];
         kk = double_to_int(pp);
         if (kk < nmin[ii])
            nmin[ii] = kk;
         if (kk > nmax[ii])
            nmax[ii] = kk;
      }
   }

   di = -nmin[0];
   dj = -nmin[1];
   dk = -nmin[2];
   sni = nmax[0] - nmin[0] + 1;
   snj = nmax[1] - nmin[1] + 1;
   snk = nmax[2] - nmin[2] + 1;
   ssfo.xx = sf.origin.xx - sf.stepv[0].xx * double(di) - sf.stepv[1].xx * double(dj) - sf.stepv[2].xx * double(dk);
   ssfo.yy = sf.origin.yy - sf.stepv[0].yy * double(di) - sf.stepv[1].yy * double(dj) - sf.stepv[2].yy * double(dk);
   ssfo.zz = sf.origin.zz - sf.stepv[0].zz * double(di) - sf.stepv[1].zz * double(dj) - sf.stepv[2].zz * double(dk);

   if (sni < 1 || snj < 1 || snk < 1)
      return -1;

   scscalar = allocate_scalar(sni, snj, snk);

   for (ii = 0; ii < sni; ii++)
   for (jj = 0; jj < snj; jj++)
   for (kk = 0; kk < snk; kk++)
   {
      if (!box_check(ii, jj, kk, &ssfo, sf.stepv, &sco, scvn, nsc))
         scscalar[ii][jj][kk] = 0.0;
      else if (ii >= di && ii <= di + sf.ni && jj >= dj && jj <= dj + sf.nj && kk >= dk && kk <= dk + sf.nk)
      {
         ui = ii - di;
         uj = jj - dj;
         uk = kk - dk;
         if (ui == sf.ni)
            ui = 0;
         if (uj == sf.nj)
            uj = 0;
         if (uk == sf.nk)
            uk = 0;
         scscalar[ii][jj][kk] = sf.scalar[ui][uj][uk];
      }
      else if (box_check(ii, jj, kk, &ssfo, sf.stepv, &uco, ucvn, nuc))
         scscalar[ii][jj][kk] = 0.0;
      else
      {
         point.xx = ssfo.xx - sf.origin.xx + sf.stepv[0].xx * double(ii) + sf.stepv[1].xx * double(jj) + sf.stepv[2].xx
	   * double(kk) - ucv[0].xx * double(nuc[0]) - ucv[1].xx * double(nuc[1]) - ucv[2].xx * double(nuc[2]);
         point.yy = ssfo.yy - sf.origin.yy + sf.stepv[0].yy * double(ii) + sf.stepv[1].yy * double(jj) + sf.stepv[2].yy
	   * double(kk) - ucv[0].yy * double(nuc[0]) - ucv[1].yy * double(nuc[1]) - ucv[2].yy * double(nuc[2]);
         point.zz = ssfo.zz - sf.origin.zz + sf.stepv[0].zz * double(ii) + sf.stepv[1].zz * double(jj) + sf.stepv[2].zz
	   * double(kk) - ucv[0].zz * double(nuc[0]) - ucv[1].zz * double(nuc[1]) - ucv[2].zz * double(nuc[2]);
         for (ll = 0; ll < 3; ll++)
         {
            pp = (point.xx * sfsn[ll].xx + point.yy * sfsn[ll].yy + point.zz * sfsn[ll].zz) / sfsw[ll];
	    nn[ll] = double_to_int(pp);
         }
         if (nn[0] >= 0 && nn[0] <= sf.ni && nn[1] >= 0 && nn[1] <= sf.nj && nn[2] >= 0 && nn[2] <= sf.nk)
         {
            ui = nn[0];
            uj = nn[1];
            uk = nn[2];
            if (ui == sf.ni)
               ui = 0;
            if (uj == sf.nj)
               uj = 0;
            if (uk == sf.nk)
               uk = 0;
            scscalar[ii][jj][kk] = sf.scalar[ui][uj][uk];
         }
         else
            scscalar[ii][jj][kk] = 0.0;
      }
   }

   if (sf.scalar != NULL)
   {
      deallocate_scalar(sf.scalar, sf.ni, sf.nj);
   }
   sf.scalar = scscalar;

   sf.ni = sni;
   sf.nj = snj;
   sf.nk = snk;

   sf.origin.xx = ssfo.xx;
   sf.origin.yy = ssfo.yy;
   sf.origin.zz = ssfo.zz;

   return 0;
}

double inversion(int nn, double *xx, double *yy, double zz)
{
   int ii, jj;
   double aa, bb, cc, ss, tt, vv, ww, uu = INF9;

   for (ii = 0; ii < nn - 2; ii++)
   {
      vv = (yy[ii] - zz) * (yy[ii] - zz) + (yy[ii + 1] - zz) * (yy[ii + 1] - zz) + (yy[ii + 2] - zz) * (yy[ii + 2] - zz);
      if (vv < uu)
      {
         jj = ii;
         uu = vv;
      }
   }

   aa = ((xx[jj + 1] - xx[jj]) * (yy[jj + 2] - yy[jj + 1]) - (xx[jj + 2] - xx[jj + 1]) * (yy[jj + 1] - yy[jj])) / ((xx[jj + 2] - xx[jj]) * (xx[jj + 2] - xx[jj + 1]) * (xx[jj + 1] - xx[jj]));
   bb = ((xx[jj + 1] - xx[jj]) * (yy[jj + 2] - yy[jj + 1]) + (xx[jj + 2] - xx[jj + 1]) * (yy[jj + 1] - yy[jj])) / (2.0 * (xx[jj + 2] - xx[jj + 1]) * (xx[jj + 1] - xx[jj])) - aa * (xx[jj] + 2.0 * xx[jj + 1] + xx[jj + 2]) / 2.0;
   cc = (yy[jj] + yy[jj + 1] + yy[jj + 2] - bb * (xx[jj] + xx[jj + 1] + xx[jj + 2]) - aa * (xx[jj] * xx[jj] + xx[jj + 1] * xx[jj + 1] + xx[jj + 2] * xx[jj + 2])) / 3.0;

   uu = (-bb + sqrt(bb * bb - 4.0 * aa * (cc - zz))) / (2.0 * aa);
   vv = (-bb - sqrt(bb * bb - 4.0 * aa * (cc - zz))) / (2.0 * aa);

   ss = (xx[jj] - uu) * (xx[jj] - uu) + (xx[jj + 1] - uu) * (xx[jj + 1] - uu) + (xx[jj + 2] - uu) * (xx[jj + 2] - uu);
   tt = (xx[jj] - vv) * (xx[jj] - vv) + (xx[jj + 1] - vv) * (xx[jj + 1] - vv) + (xx[jj + 2] - vv) * (xx[jj + 2] - vv);

   if (ss < tt)
      ww = uu;
   else
      ww = vv;

   return ww;
}

double isovalue_scale(int power, double isovalue, const scalar_field & sff)
{
   const int ns = 64;
   int ii, jj, kk, pp, ss;
   double samp, smin = INF9, smax = -INF9;
   double sa[ns], sb[ns], sc[ns], sd, se, sf;

   for (ii = 0; ii < sff.ni; ii++)
   for (jj = 0; jj < sff.nj; jj++)
   for (kk = 0; kk < sff.nk; kk++)
   {
      if (sff.scalar[ii][jj][kk] < smin)
         smin = sff.scalar[ii][jj][kk];
      if (sff.scalar[ii][jj][kk] > smax)
         smax = sff.scalar[ii][jj][kk];
   }
   if (fabs(smax) > fabs(smin))
      samp = fabs(smax);
   else
      samp = fabs(smin);

   if (power > 0)
   {
      for (ss = 0; ss < ns; ss++)
         sa[ss] = double(ss) / double(ns - 1);
      for (ss = 0; ss < ns; ss++)
         sb[ss] = sa[ss] * samp;
      sb[0] = sb[0] - EPS9;
      sb[ns - 1] = sb[ns - 1] + EPS9;
      for (ss = 0; ss < ns; ss++)
         sc[ss] = 0.0;
      for (ii = 0; ii < sff.ni; ii++)
      for (jj = 0; jj < sff.nj; jj++)
      for (kk = 0; kk < sff.nk; kk++)
      {
         sd = fabs(sff.scalar[ii][jj][kk]);
         se = 1.0;
         for (pp = 0; pp < power; pp++)
            se = se * sd;
         for (ss = 0; ss < ns; ss++)
            if (sb[ss] < sd)
               sc[ss] = sc[ss] + se;
      }
      sf = sc[0];
      if (sf > EPS9)
      {
         for (ss = 0; ss < ns; ss++)
            sc[ss] = sc[ss] / sf;
         isovalue = inversion(ns, sa, sc, isovalue);
      }
      else
         isovalue = 0.0;
   }
   isovalue = isovalue * samp;

   return isovalue;
}

// round to nearest integer
int double_to_int(double dd)
{
   if (dd < 0.0)
      return int(dd - 0.5);
   else
      return int(dd + 0.5);
}

