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
#include "inc.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>



void options::defaults() {
 mouse    = false;
 search   = string("");
 protein  = string("");
 date     = string("");
 dIRZ            = 0.1;
 dIRO            = 1;
 regions[0]      = 1000;
 regions[1]      = 2000;
 regions[2]      = 3000;
 regions[3]      = 4000;
 inExon          = 50;
 inIntron        = 200;
 mlen            = 0;
 microarray_rows = 0;
 path_splicing_change = string("");
 tetramers       = string("");
}

options::options() {
 defaults();
}

void options::print(ostream & out) {
 cout << "Mouse:        " << mouse << endl;
 cout << "Analysis folder:      " << protein<< endl;
 cout << "Search RSYW:  " << search << endl;
 cout << "dIRZ:         " << dIRZ << endl;
 cout << "dIRO:         " << dIRO << endl;
 cout << "inExon:       " << inExon << endl;
 cout << "inIntron:     " << inIntron << endl;
 cout << "Path SP:      " << path_splicing_change << endl;
 cout << "Path Tetra:   " << tetramers << endl;
 cout << "Array rows:   " << microarray_rows << endl << endl;
 cout << "Motif Length: " << mlen << endl;
}

options::options(int argc,char *argv[]) {
 defaults();
 AnyOption opt;
 opt.addUsage( "tetramer " );
 opt.addUsage( "Usage: " );
 opt.addUsage( "" );
 opt.addUsage( " -h  --help             Prints this help " );
 opt.addUsage( " -c  --config-file      Read configuration file" );
	
 opt.addUsage( "" );
 opt.setOption( "config-file", 'c');

 opt.setFlag(  "help", 'h' );

 opt.processCommandArgs( argc, argv );

 if ( ! opt.hasOptions()) opt.printUsage();
 if ( opt.getFlag( "help" ) || opt.getFlag( 'h' ) ) opt.printUsage();
 if( opt.getValue( 'c' ) != NULL  || opt.getValue( "config-file" ) != NULL  ){
  
  	ifstream in(opt.getValue( 'c' ));
  	
  	  if(in.is_open()){
  	  		string line;
  	  		int i = 0;
   			while(getline(in,line)){
   				cout << line.c_str() << endl;
   				string param = line.substr(line.find("=")+1,line.length());
				cout << param.c_str() << endl;
				if(i==0){
					if(strcmp(param.c_str(),"Mouse")==0) mouse=true;
				} else if (i==1){
					path_splicing_change = param;
				} else if (i==2){
					tetramers = param;
				} else if (i==3){
					protein=param;
				} else if (i==4){
					search=param;
				}
				
   				i++;
   			}
			  if (search.compare("R")==0){
					tetramers = tetramers + string("/r/");	
					protein   = protein + string("/r/");
			  }else if (search.compare("N")==0) {
				   tetramers = tetramers + string("/nr/");
				   protein   = protein + string("/nr/");
			}

			 int ret = mkdir(protein.c_str(),0775);
			 microarray_rows = (mouse)?(32874):(53624);
  			
  		}else {
  			cout << "invalid configuration file" << endl;
  		}
  }
}

//------------------------------------------------------

bln::bln(){}

bln::bln(string &line){
 istringstream linestream(line);
 string token;
 vector<string> row;
 while (getline(linestream, token, '\t')) row.push_back(token);
 if(row.size()>0){
  strcpy(chr,row[0].c_str());
  chromStart = atoi(row[1].c_str());
  chromEnd   = atoi(row[2].c_str());
  strand_num = atoi(row[3].c_str());
  forw = (strand_num>0)?(true):(false);
  if(pos_map_end && cat_of_change) sprintf(mix,"CE_%d_%d",cat_of_change,pos_map_end);
  cDNA = 1;
 }
}

bln::bln(const bln & tocopy){
 rowID         = tocopy.rowID;
 myRID         = tocopy.myRID;
 cat_of_change = tocopy.cat_of_change;
 strcpy(chr,tocopy.chr);
 strand_num    = tocopy.strand_num;
 forw          = tocopy.forw;
 cDNA          = tocopy.cDNA;
 chromStart    = tocopy.chromStart;
 chromEnd      = tocopy.chromEnd;
 pos_map_start = tocopy.pos_map_start;
 pos_map_end   = tocopy.pos_map_end;
 strcpy(mix,tocopy.mix);
}

void bln::print(ostream &out){
 out << chr << '\t'
     << chromEnd << '\t'
     << pos_map_end << '\t'
     << strand_num << '\t'
     << myRID  << '\t'
     << rowID  << '\t'
     << cat_of_change << endl;
};

void bln::printAll(ostream &out){
 out << chr << '\t'
     << chromStart << '\t'
     << chromEnd << '\t'
     << strand_num << '\t';
 if(forw)
  out << pos_map_start << '\t'
      << pos_map_end << '\t';
 else
  out << pos_map_end << '\t'
      << pos_map_start << '\t';
 out << myRID  << '\t'
     << rowID  << '\t'
     << cat_of_change << endl;
};

void bln::printRank(ostream &out){
 out << pos_map_end << '\t'
     << cat_of_change << '\t'
     << cDNA << endl;
};

void bln::readExpandedLine(string &line){}

void bln::set_bed_on_the_map ( options opt, bool forw, int region, unsigned int region_start, unsigned int region_stop,
                               unsigned int geno_bound, unsigned int inIntron, unsigned int inExon ){

 int delta_start = abs( (int)geno_bound - (int)chromStart ),
     delta_stop  = abs( (int)chromEnd - (int)geno_bound );

 // [ region_start, region_stop ]
 if( region_start <= chromStart && chromEnd <= region_stop ){
  if(forw){
   if( chromStart <= geno_bound) pos_map_start = opt.regions[ region ] -  delta_start ;
   else                          pos_map_start = opt.regions[ region ] +  delta_start ;

   if( chromEnd <= geno_bound)   pos_map_end   = opt.regions[ region ] -  delta_stop ;
   else                          pos_map_end   = opt.regions[ region ] +  delta_stop ;
  }else{
   if(chromStart <= geno_bound)  pos_map_start = opt.regions[ 3-region ] + delta_start;
   else                          pos_map_start = opt.regions[ 3-region ] - delta_start;

   if(chromEnd <= geno_bound)    pos_map_end   = opt.regions[ 3-region ] + delta_stop;
   else                          pos_map_end   = opt.regions[ 3-region ] - delta_stop;
  }

  // [region_start, no]
  }else if(region_start<=chromStart && chromEnd > region_stop ){
   if(forw){
    if(chromStart <= geno_bound) pos_map_start = opt.regions[ region ] -  delta_start ;
    else                         pos_map_start = opt.regions[ region ] +  delta_start ;

    if(region==0 || region==2)   pos_map_end = opt.regions[ region ] + inIntron ;  // end
    else                         pos_map_end = opt.regions[ region ] + inExon ;  // end
   }else{

    if(chromStart <= geno_bound) pos_map_start = opt.regions[ 3-region ] + delta_start;
    else                         pos_map_start = opt.regions[ 3-region ] - delta_start;

    if(region==0 || region==2)   pos_map_end = opt.regions[ 3-region ] - inIntron ;  // end
    else                         pos_map_end = opt.regions[ 3-region ] - inExon ;  // end

   }
   // [ no, firstaregion_stop]
   }else if(chromStart<region_start && chromEnd<=region_stop){
    if(forw){
     if(region==0 || region==2)  pos_map_start = opt.regions[ region ] - inExon;  // start
     else                        pos_map_start = opt.regions[ region ] - inIntron;

     if(chromEnd <= geno_bound)  pos_map_end   = opt.regions[ region ] -  delta_stop ;
     else                        pos_map_end   = opt.regions[ region ] +  delta_stop ;
    }else{

     if(region==0 || region==2)  pos_map_start = opt.regions[ 3-region ] + inExon;
     else                        pos_map_start = opt.regions[ 3-region ] + inIntron;

     if(chromEnd <= geno_bound) pos_map_end   = opt.regions[ 3-region ] + delta_stop;
     else                       pos_map_end   = opt.regions[ 3-region ] - delta_stop;
    }

   // [ no, no ]
   }else if( chromStart<region_start && chromEnd>region_stop){
    if(forw){
     if(region==0 || region==2){
      pos_map_start = opt.regions[ region ] - inExon;
      pos_map_end   = opt.regions[ region ] + inIntron;
     }else{
      pos_map_start = opt.regions[ region ] - inIntron;
      pos_map_end   = opt.regions[ region ] + inExon;
     }
    }else{
    if(region==0 || region==2){
      pos_map_start = opt.regions[ 3-region ] + inExon;
      pos_map_end   = opt.regions[ 3-region ] - inIntron;
     }else{
      pos_map_start = opt.regions[ 3-region ] + inIntron;
      pos_map_end   = opt.regions[ 3-region ] - inExon;
     }
    }
   }
};


void bln::set_bed_middle_on_the_map( options opt, int middle_position, bool forw, int region, unsigned int region_start, unsigned int region_stop, unsigned int geno_bound, unsigned int inIntron, unsigned int inExon ){

 int delta = abs( (int)geno_bound - (middle_position) );

 int pos_map;
 if(forw){
  if( (unsigned int) middle_position <= geno_bound) pos_map = opt.regions[ region ] -  delta ;
  else                               pos_map = opt.regions[ region ] +  delta ;
 }else{
  if((unsigned int) middle_position <= geno_bound)  pos_map = opt.regions[ 3-region ] + delta;
  else                               pos_map = opt.regions[ 3-region ] - delta;
 }
 pos_map_start = pos_map;
 pos_map_end   = pos_map;
};

bool sortingBed(const bln & b1,const bln & b2){
 return b1.rowID<b2.rowID;
}

bool sortingMyRID(const bln & b1,const bln & b2){
 return b1.myRID<b2.myRID;
}

bool sortingCOC(const bln & b1,const bln & b2){
 return b1.cat_of_change<b2.cat_of_change;
}

bool sortingPosMap(const bln & b1,const bln & b2){
 return b1.pos_map_end<b2.pos_map_end;
}

bool sortingBedMix(const bln & b1,const bln & b2){
//  return (string(b1.mix) < string(b2.mix));
 int pr1=0,pr2=0;
 if(b1.cat_of_change==1)         pr1=10;
 else if(b1.cat_of_change==0)    pr1=1;
 else if(b1.cat_of_change==(-1)) pr1=(-10);
 if(b2.cat_of_change==1)         pr2=10;
 else if(b2.cat_of_change==0)    pr2=1;
 else if(b2.cat_of_change==(-1)) pr2=(-10);
 int a = pr1*b1.pos_map_end;
 int b = pr2*b2.pos_map_end;
 if(a<0 && b<0) return a>b;
 else  return a<b; //quasi ok
}
//------------------------------------------------------

windows::windows(){};

windows::windows(const bln &tocopy):bln(tocopy){
 sum1=0;
 sum2=0;
 sum3=0;
 sum4=0;
}

void windows::setSums(bln mybln,options opt){
 unsigned int reg1_start = opt.regions[0]-opt.inExon,   reg1_stop  = opt.regions[0]+opt.inIntron,
              reg2_start = opt.regions[1]-opt.inIntron, reg2_stop  = opt.regions[1]+opt.inExon,
              reg3_start = opt.regions[2]-opt.inExon,   reg3_stop  = opt.regions[2]+opt.inIntron,
              reg4_start = opt.regions[3]-opt.inIntron, reg4_stop  = opt.regions[3]+opt.inExon;

 if(reg1_start<=mybln.pos_map_end && mybln.pos_map_end<=reg1_stop)      sum1+=mybln.cDNA;
 else if(reg2_start<=mybln.pos_map_end && mybln.pos_map_end<=reg2_stop) sum2+=mybln.cDNA;
 else if(reg3_start<=mybln.pos_map_end && mybln.pos_map_end<=reg3_stop) sum3+=mybln.cDNA;
 else if(reg4_start<=mybln.pos_map_end && mybln.pos_map_end<=reg4_stop) sum4+=mybln.cDNA;
}

void windows::print_myRID(ofstream &out){
 out << myRID << "\t" << chromStart << "\t" << sum1 << "\t" << sum2 << "\t" << sum3 << "\t" << sum4 << endl;
}

void windows::print_rowID(ofstream &out){
 out << rowID << "\t" << chromStart << "\t" << sum1 << "\t" << sum2 << "\t" << sum3 << "\t" << sum4 << endl;
}

//------------------------------------------------------

void expand(vector<bln> bedres, vector<bln> &bedexp){
 if(bedres.size()>0){
  unsigned int gen_end,rna_end;
  for(unsigned int g=0; g<bedres.size();g++){
   gen_end = bedres[g].chromEnd;
   rna_end = bedres[g].pos_map_end;
   if(bedres[g].strand_num>0){
    for(int i=0; i<bedres[g].strand_num; i++){
     bln bb(bedres[g]);
     bb.chromEnd    = gen_end - i;
     bb.pos_map_end = rna_end - i;
     sprintf(bb.mix,"CE_%d_%d",bb.pos_map_end,bb.cat_of_change);
     bedexp.push_back(bb);
    }
   } else {
    for( int i=0; i<(-bedres[g].strand_num); i++){
     bln bb(bedres[g]);
     bb.chromEnd    = gen_end - i;
     bb.pos_map_end = rna_end + i;
     sprintf(bb.mix,"CE_%d_%d",bb.pos_map_end,bb.cat_of_change);
     bedexp.push_back(bb);
    }
   }
  }
 }else{
  cout << "no size" << endl;
  throw ("Expand: no results");
 }
}

void dIrank( vector<bln> bedexp, options opt, vector <windows> &wins){
 unsigned int c=0,a=0;
 wins.push_back( windows(bedexp[a]) );
 wins[c].setSums( bedexp[a], opt );
 for( a=1; a<bedexp.size(); a++ ){
  if( bedexp[a].myRID == bedexp[a-1].myRID ){
   wins[c].setSums(bedexp[a],opt);
  }else{
   c++;
   wins.push_back( windows( bedexp[a] ) );
   wins[c].setSums( bedexp[a],opt );
  }
 }
}

//------------------------------------------------------
int EM_clip_microarray_CE(options &opt,const char *sp_fname,const char *bed_fname, const char *output_fname, const char * tet_name,ofstream & stats){
 try{
  ifstream in(sp_fname), fbed(bed_fname);
  ofstream fout(output_fname),fstats(),fres1("/tmp/out.txt");
  string line,f;
  vector<bln> bedlines, bedres;
  unsigned int CEone      = 0, rCEone    = 0,
               CEminone   = 0, rCEminone = 0,
               CEzero     = 0, rCEzero   = 0,
               offset     = 0;

  if(in.is_open() && fbed.is_open()){
   while(getline(fbed,line)) bedlines.push_back(bln(line)); // Read BED file
   fbed.close();

   while(getline(in,f)){  // read splicing change lines and evaluate bed positions

    gString spline=f;
    vector<gString> pars;
    spline.split(';',pars);

    string       chr        = string(pars[2]);
    bool         forw       = (string(pars[3]).compare("+")==0)?(true):(false);
    unsigned int myRID      = atoi(string(pars[0]).c_str()),
                 rowID      = atoi(string(pars[1]).c_str()),
                 skip_start = atoi(string(pars[4]).c_str()) - offset, // posizioni genomiche
                 in_start   = atoi(string(pars[5]).c_str()),
                 in_stop    = atoi(string(pars[6]).c_str()),
                 skip_stop  = atoi(string(pars[7]).c_str()) + offset;
    double       dIRank     = atof(string(pars[8]).c_str());
    int          cofc       = 2;

    unsigned int firstaregion_start  = skip_start - opt.inExon,
                 firstaregion_stop   = skip_start + opt.inIntron,
                 secondaregion_start = in_start - opt.inIntron,
                 secondaregion_stop  = in_start + opt.inExon,
                 thirdaregion_start  = in_stop - opt.inExon,
                 thirdaregion_stop   = in_stop + opt.inIntron,
                 fourtharegion_start = skip_stop - opt.inIntron,
                 fourtharegion_stop  = skip_stop + opt.inExon,
                 left_intron_middle  = skip_start + round((in_start-skip_start)/2),
                 exon_middle         = in_start   + round((in_stop-in_start)/2),
                 right_intron_middle = in_stop    + round((skip_stop-in_stop)/2);

    if((-opt.dIRZ)<=dIRank && dIRank<=opt.dIRZ){
     CEzero++;
     cofc = 0;
    }else if(dIRank>=opt.dIRO){
     CEone++;
     cofc = 1;
    }else if(dIRank<=(-opt.dIRO)){
     CEminone++;
     cofc= -1;
    }
    /*----------------------------
      updating regions boundaries
    -----------------------------*/
    if( left_intron_middle < firstaregion_stop )   firstaregion_stop   = left_intron_middle;
    if( left_intron_middle > secondaregion_start ) secondaregion_start = left_intron_middle;

    if( exon_middle < secondaregion_stop ) secondaregion_stop = exon_middle;
    if( exon_middle > thirdaregion_start ) thirdaregion_start = exon_middle;

    if( right_intron_middle < thirdaregion_stop )   thirdaregion_stop   = right_intron_middle;
    if( right_intron_middle > fourtharegion_start ) fourtharegion_start = right_intron_middle;

    unsigned int inIntronLeft  = firstaregion_stop - skip_start,
                 inExon        = secondaregion_stop - in_start,
                 inIntronRight = thirdaregion_stop - in_stop;
    bool spara=false;
    for(unsigned int g=0; g<bedlines.size();g++){
     if(strcmp(chr.c_str(),bedlines[g].chr)==0 &&
       (forw==bedlines[g].forw) && 
       (bedlines[g].chromStart<=fourtharegion_stop) && (bedlines[g].chromEnd>=firstaregion_start) ){
      bln bed(bedlines[g]);
      bed.rowID         = rowID;
      bed.myRID         = myRID;
      bed.cat_of_change = cofc;
      /*--------------------
           FIRST REGION
      --------------------*/
      if( bedlines[g].chromStart <= firstaregion_stop && bedlines[g].chromEnd >= firstaregion_start ){
       bed.set_bed_on_the_map(opt, forw, 0, firstaregion_start, firstaregion_stop, skip_start, inIntronLeft, opt.inExon);
       bedres.push_back(bed);
       //spara=true;
      }
      /*--------------------
           SECOND REGION
      --------------------*/
      else if( bedlines[g].chromStart <= secondaregion_stop && bedlines[g].chromEnd >= secondaregion_start ){
       bed.set_bed_on_the_map(opt, forw, 1, secondaregion_start, secondaregion_stop, in_start, inIntronLeft, inExon);
       bedres.push_back(bed);
       if(bedlines[g].chromEnd >= (secondaregion_start+100) && bedlines[g].chromStart <= (secondaregion_stop-15)) spara=true;
      }
      /*--------------------
          THIRD REGION
      --------------------*/
      else if( bedlines[g].chromStart <= thirdaregion_stop && bedlines[g].chromEnd >= thirdaregion_start ){
       bed.set_bed_on_the_map(opt, forw, 2, thirdaregion_start, thirdaregion_stop, in_stop, inIntronRight, inExon);
       bedres.push_back(bed);
       if(bedlines[g].chromEnd >= (thirdaregion_start+15) && bedlines[g].chromStart <= (thirdaregion_stop-130)) spara=true;
      }
      /*--------------------
           FOURTH REGION
      --------------------*/
      else if( bedlines[g].chromStart <= fourtharegion_stop && bedlines[g].chromEnd >= fourtharegion_start ){
       bed.set_bed_on_the_map(opt, forw, 3, fourtharegion_start, fourtharegion_stop, skip_stop, inIntronRight, opt.inExon);
       bedres.push_back(bed);
       //spara=true;
      }
     }
    }
    if(spara==true){
     if(cofc==0){          rCEzero++;
     }else if(cofc==1){    rCEone++;
     }else if(cofc==(-1)){ rCEminone++;
     }
    }
  } // end splicing change file
  fout << "CE1 = " << CEone << ", CE-1 = " << CEminone <<", CE0 = " << CEzero << endl;
  stats << tet_name << "\t" << rCEone << "\t" << rCEminone << "\t" << rCEzero << endl;

  if(bedres.size()>0){

   // Expand
   // -------
   vector<bln> bedexp;
   expand( bedres, bedexp );
   sort(bedexp.begin(),bedexp.end(),sortingBedMix); //    for(unsigned int g=0; g<bedexp.size();g++) bedexp[g].print(fres1);
   
   
   // dIrank
   //-------
   vector <bln> bedrank;
   unsigned int c=0, a=0;
   while( a<bedexp.size() ){
    bln bb(bedexp[a]);
    if(a==0){
     bb.cDNA=1;
     bedrank.push_back( bb );
    }else{
     if(bedexp[a].pos_map_end==bedexp[a-1].pos_map_end && bedexp[a].cat_of_change==bedexp[a-1].cat_of_change){
      bedrank[c].cDNA++;
     }else{
      c++;
      bb.cDNA=1;
      bedrank.push_back( bb );
     }
    }
    a++;
   }

   for(unsigned int g=0; g<bedrank.size();g++) bedrank[g].printRank(fout);

   }else cout << "no results!" << endl;
  }else cout << "Impossible open data file: "<< sp_fname << endl;
 }catch(gException &e){
  cout << e.what() << endl;
 }catch(...){ 
  cout <<"\tException";
 }
 return 0;
}

int count_per_regions(options &opt,const char *sp_fname,const char *bed_fname, const char *output_fname){
	try{
		ifstream in(sp_fname), fbed(bed_fname);
		ofstream fout(output_fname);
		string line,f;
		vector<bln> bedlines, bedres;
		unsigned int CEone      = 0, 
					 CEminone   = 0, 
    	  			 CEzero     = 0;
		
		fout << "myRID\trowID\ttype\thits_region1\thits_region2\thits_region3\n"; 
		if(in.is_open() && fbed.is_open()){
			
			while(getline(fbed,line)) bedlines.push_back(bln(line)); // Read BED file
			fbed.close();
			
			while(getline(in,f)){  // read splicing change lines and evaluate bed positions

				/*----------------------------
				  LOADING EXON
				 -----------------------------*/
				gString spline=f; 

				vector<gString> pars;
				spline.split(';',pars);
				string       chr        = string(pars[2]);
				bool         forw       = (string(pars[3]).compare("+")==0)?(true):(false);
				unsigned int myRID      = atoi(string(pars[0]).c_str()),
							 rowID      = atoi(string(pars[1]).c_str()),
							 in_start   = atoi(string(pars[5]).c_str()),
							 in_stop    = atoi(string(pars[6]).c_str());
				double       dIRank     = atof(string(pars[8]).c_str());
				
				unsigned int exon_region1_stop   = in_start + 30,
							 exon_region2_start  = in_stop  - 30,
							 exon_middle         = in_start + round((in_stop-in_start)/2), 
							 r1_start            = in_start - 35,
							 r1_stop             = in_start - 5,
							 r3_start            = in_stop  + 10,
							 r3_stop             = in_stop  + 40;

				
				int  cofc       = 2;
				if((-opt.dIRZ)<=dIRank && dIRank<=opt.dIRZ){
					CEzero++;
					cofc = 0;
				}else if(dIRank>=opt.dIRO){
					CEone++;
					cofc = 1;
				}else if(dIRank<=(-opt.dIRO)){
					CEminone++;
					cofc= -1;
				}
				
				/*----------------------------
				 updating regions boundaries
				 -----------------------------*/
			
				if( exon_middle < exon_region1_stop )  exon_region1_stop  = exon_middle;
				if( exon_middle > exon_region2_start ) exon_region2_start = exon_middle;
			

				bool hits_region1 = false,
					 hits_region2 = false,
					 hits_region3 = false;

				/*----------------------------
				 load teramer hits 
				 -----------------------------*/
				
				for(unsigned int g=0; g<bedlines.size();g++){
					if(strcmp(chr.c_str(),bedlines[g].chr)==0 &&
					   (forw==bedlines[g].forw)){ 

							if (!forw){
								
								
								r1_start            = in_stop + 5, 
								r1_stop             = in_stop + 35, // piu alta
								r3_start            = in_start - 40, // piu bassa
								r3_stop             = in_start - 10; 
								
								
								
							}
							
							// REGION 1 
							if( bedlines[g].chromStart <= r1_stop && bedlines[g].chromEnd >= r1_start )	{ hits_region1 = true; }
							// REGION 2 
							else if( bedlines[g].chromStart <= exon_region1_stop && bedlines[g].chromEnd >= in_start ){ hits_region2 = true; }
							else if( bedlines[g].chromStart <= in_stop && bedlines[g].chromEnd >=exon_region2_start  ){	hits_region2 = true; }
							// REGION 3
							else if( bedlines[g].chromStart <= r3_stop && bedlines[g].chromEnd >= r3_start ){ hits_region3 = true; }
				
					}
				}

				fout << myRID << "\t" << rowID  << "\t" << cofc  << "\t"  << hits_region1  << "\t" <<  hits_region2 << "\t" << hits_region3 << "\n"; 
				
			}
			
		}else cout << "Impossible open data file: "<< sp_fname << endl;
	}catch(gException &e){
		cout << e.what() << endl;
	}catch(...){ 
		cout <<"\tException";
	}
	return 0;
}

int counting_per_region(options &opt){
	
	string folder  = opt.protein,		
		   command = string("wc -l ") + opt.tetramers + string("*.bed > ") + opt.tetramers + string("tetramers_files.txt");
	
	system(command.c_str());
	
	string fin = opt.tetramers + string("tetramers_files.txt");
	
	cout << fin.c_str() << endl;
	
	ifstream in(fin.c_str());
	
	if(in.is_open()){ // load tetramers bed file in the folder
		
		unsigned int rows_bed;
		string bed_fname, output_fname,
		splicing_fname = opt.path_splicing_change,
		filelist_name  = folder + string("filelist_count.tsv");
		
		ofstream filelist(filelist_name.c_str());
		
		while(!in.eof()){ // load single tet bed file
			
			in >> rows_bed >> bed_fname;
			if(bed_fname.compare("total")!=0){
				size_t p1       = bed_fname.find_last_of('/');
				string tet_name = bed_fname.substr(p1+1).substr(0,4);
				if(bed_fname.substr(p1+1).compare("totale")!=0){
					
					string fout =  opt.protein + tet_name.c_str() + string("_region_count.tsv");
					
					count_per_regions(opt,splicing_fname.c_str(),bed_fname.c_str(),fout.c_str());
					
					filelist << tet_name.c_str() << string("_region_count.tsv") << endl;
					
				}
			}
		}
	}else cout << "impossible opening file" << endl;
	in.close();
	
	//loop on splicing change human files
	return 0;
}

int tetramer(options &opt){	 
	/*-------------------------------
	 
	 RUN ANALISYS
	 
	 -------------------------------*/
	int ret = 0;
	
	string folder  =  opt.protein ,
		command = string("wc -l ") + opt.tetramers + string("*.bed > ") + opt.tetramers + string("tetramers_files.txt");

		system(command.c_str());

		string fin =  opt.tetramers + string("tetramers_files.txt");
		ifstream in(fin.c_str());
		
		if(in.is_open()){
			string bed_fname, output_fname,
			splicing_fname = opt.path_splicing_change ,
			stats_name     = folder + string("STATS.txt"),
			filelist_name  = folder + string("filelist.txt");
			
			unsigned int rows_bed;
			cout << stats_name << endl;
			ofstream stats(stats_name.c_str()),
			filelist(filelist_name.c_str());
			while(!in.eof()){
				in >> rows_bed >> bed_fname;
				if(bed_fname.compare("total")!=0){
					size_t p1       = bed_fname.find_last_of('/');
					string tet_name = bed_fname.substr(p1+1).substr(0,4);
					if(bed_fname.substr(p1+1).compare("totale")!=0){
						
						string fout = opt.protein + bed_fname.substr(p1+1).c_str();
						
						ret = EM_clip_microarray_CE(opt,splicing_fname.c_str(),bed_fname.c_str(),fout.c_str(),tet_name.c_str(),stats);
						
						filelist << bed_fname.substr(p1+1).c_str() << endl;
						
					}
				}
			}
		}else cout << "impossible opening file" << endl;
		in.close();
	return ret;
}

