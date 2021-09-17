#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

typedef struct {
   int number;
   char symbol[4];
   char name[14];
   int nvelec;
   int valence;
   int period;
   int group;
   double rcov;
   double rvdw;
   double mass;
   double color[3];
} ELEMENT;

int get_num_elements();
void get_symbol_by_z(int *zz, char symbol[3]);
void get_symbol_by_z_(int *zz, char symbol[3]);

#endif

