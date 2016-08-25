/***************************************************************************
 *   Copyright (C) 2010 by Matteo Cereda   *
 *   mcereda@tinca.bp.lnf.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "struct.h"
#include "anyoption.h"

#include "geco_define.h"
#include "geco_base.h"
#include "geco_array.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <typeinfo>

#include <vector>
#include <exception>

using namespace std;
using namespace geco;


// Common Functions
//-----------------------
int splita(char *original, char * split1, char * split2);
int splitb(char *original, char * split1, char * split2);
int compare_ints( int a, int b );


class options {
public:
    bool         mouse;
    string       search;
    string       protein;
    string       date;
    string       date_bed;
    double       dIRZ;
    double       dIRO;
    unsigned int regions[4];
    unsigned int inExon;
    unsigned int inIntron;
    string       path;
    string       path_splicing_change;
    string       tetramers;
    string       path_fbed;
    string       dir;
    string       percent;
    double        mlen;
    unsigned int microarray_rows; // mouse 32874 human 53624
    options();
    options(int argc,char *argv[]);
    void defaults();
    void print(ostream & out);
};

class search_opt {
public:
    bool         mouse;
    double        mlen;
    string       date;
    string       date_bed;
    string       path;
    string       tetramers;
    string       search;
    string       path_fbed;
    search_opt();
    search_opt(int argc,char *argv[]);
    void defaults();
    void print(ostream & out);
};

class bln{
 public:
  char          chr[10];
  unsigned int  chromStart;
  unsigned int  chromEnd;
  unsigned int  pos_map_start;
  unsigned int  pos_map_end;
  int           strand_num;
  bool          forw;
  unsigned int  rowID;
  unsigned int  myRID; // relative rowID to CE splicing files.
  int           cat_of_change;
  char          mix[40];
  int           cDNA;

  bln();
  bln(const bln &tocopy);
  bln(string &line);

  bool operator <(const bln &other) const{ return rowID<other.rowID;};
  void readExpandedLine(string &line);
  void print(ostream &out);
  void printAll(ostream &out);
  void printRank(ostream &out);

  void set_bed_on_the_map ( options opt, bool forw, int region,
                            unsigned int region_start, unsigned int region_end,
                            unsigned int geno_bound, unsigned int inIntron, unsigned int inExon);
  void set_bed_middle_on_the_map( options opt, int middle_position, bool forw, int region, 
                            unsigned int region_start, unsigned int region_stop, 
                            unsigned int geno_bound, unsigned int inIntron, unsigned int inExon );
};


class windows:public bln{
 public:
  unsigned int sum1;
  unsigned int sum2;
  unsigned int sum3;
  unsigned int sum4;

  windows();
  windows(const bln &tocopy);

  void setSums(bln mybln,options opt);
  void print_myRID(ofstream &out);
  void print_rowID(ofstream &out);
};

bool sortingBed(const bln & b1,const bln & b2);
bool sortingMyRID(const bln & b1,const bln & b2);
bool sortingCOC(const bln & b1,const bln & b2);
bool sortingPosMap(const bln & b1,const bln & b2);
bool sortingBedMix(const bln & b1,const bln & b2);

void expand(vector<bln> bedres, vector<bln> &bedexp);
void dIrank(vector<bln> bedexp, options opt, vector <windows> &wins);

int  EM_clip_microarray_CE(options &opt,const char *sp_fname,const char *bed_fname, const char *output_fname, const char * tet_name, ofstream &stats);

int tetramer(options &opt, string &splicing_change, string &percent);

int counting_per_region(options &opt, string &splicing_change, string &percent);

int count_per_regions(options &opt,const char *sp_fname,const char *bed_fname, const char *output_fname);

int get_exon_coordinates(options &opt,const char *sp_fname, const char *output_fname);

void print_exon_coordindates(options &opt, string &splicing_fname, string &fout);

template <class T> ostream & operator << (ostream & out, const gArray<T> & array) {
 size_t lfield;
 if(typeid(T) == typeid(gScore) ){
  lfield=7;
  cout.precision(1);
  cout << fixed;
 }else{
  lfield=5;
 }
 string sep(4+lfield*array.getSize(),'-');
 out << "Pos:";
 for (gPos i=0;i<array.getSize();i++) {
  out << setw(lfield) << i;
 }
 out << endl << sep.c_str() << endl;
 out << "Val:";
 for (gPos i=0;i<array.getSize();i++) {
  if(array.isNA(i)) out << setw(lfield) << "*";
  else out << setw(lfield) << array[i];
 }
 out << endl;
 return out;
}

