// struct output3{
//  char tc[40];
//  char eex[40];
//  char ein[40];
//  char strand[2];
//  char ensemble_gene_id[40];
//  char exon_num_eex[40];
//  int exon_length;
//  char exon_num_ein[40];
//  char gene_symbol[40];
//  char in[100];
//  char skip[100];
//  char stringency[100];
//  char chr[10];
//  int rowID;
//  int in_skip_number;
//  char type[10];
//  char one[40];
//  char two[40];
//  char three[40];
//  char four[40];
//  char ones[40];
//  char twos[40];
//  char threes[40];
//  char fours[40];
//  int in_start;
//  int in_stop;
//  int skip_start;
//  int skip_stop;
//  int right_intron_middle;
//  int left_intron_middle;
//  int exon_middle;
//  int firstaregion_start;
//  int firstaregion_stop;
//  int secondaregion_start;
//  int secondaregion_stop;
//  int thirdaregion_start;
//  int thirdaregion_stop;
//  int fourtharegion_start;
//  int fourtharegion_stop;
//  int first_reg;
//  int second_reg;
//  int third_reg;
//  int fourth_reg;
//  double dIRank;
//  char dIRankstr[40];
// };


struct output3a{
 char strand[2];
 char chr[10];
 int position;
 int position_1;
 double strand_num;
 int post_rna_map;
 int post_rna_map_1;
 double norm1;
 double norm2;
 char type[10];
 int rowID;
 int cat_of_change;
 int check;
 double max;
 double FDR;
 char FDRstr[40];
 struct output3a *next;
};

struct output3gh{
 char strand[2];
 char chr[10];
 int position;
 double strand_num;
 int position_1;
 int post_rna_map;
 int post_rna_map_1;
 double norm1;
 double norm2;
 char type[10];
 int rowID;
 int cat_of_change;
 int check;
 double max;
 double FDR;
 char FDRstr[40];
};

struct output3w{
 char strand[2];
 char chr[10];
 int position;
 double strand_num;
 int position_1;
 int post_rna_map;
 int post_rna_map_1;
 double norm1;
 double norm2;
 char type[10];
 int rowID;
 int cat_of_change;
 int check;
 double max;
 double FDR;
 char FDRstr[40]; 
 };

struct output3b{
 char strand[2];
 char chr[10];
 int position;
 double strand_num;
 int post_rna_map;
 double norm1;
 double norm2;
 char type[10];
 int rowID;
 int cat_of_change;
 int check;
 double max;
 double FDR;
 char FDRstr[40];
};

struct a_dIRank{
 char strand[2];
 double strand_num;
 int post_rna_map;
 double norm1;
 double norm2;
 char type[10];
 int cat_of_change;
 char mix[40];
 int cDNA;
};

struct b_dIRank{
 char strand[2];
 double strand_num;
 int post_rna_map;
 double norm1;
 double norm2;
 char type[10];
 int cat_of_change;
 int cDNA;
};

struct output3{
 char tc[40];
 char eex[40];
 char ein[40];
 char strand[2];
 char ensemble_gene_id[40];
 char exon_num_eex[40];
 int exon_length;
 char exon_num_ein[40];
 double hits;
 char description[5000];
 char gene_symbol[40];
 char reminder[40];
 char in[100];
 char skip[100];
 char stringency[100];
 char chr[10];
 char chr_hg18[10];
 char chr_hg19[10];
 int strand_hg18;
 int strand_hg19;
 int rowID;
 int in_start_hg18;
 int in_stop_hg18;
 int skip_start_hg18;
 int skip_stop_hg18;
 int in_start_hg19;
 int in_stop_hg19;
 int skip_start_hg19;
 int skip_stop_hg19;
 int in_skip_number;
 int num1;
 int num2;
 int num3;
 char A5SS[2];
 char A3SS[2];
 char CE[2];
 char MXE[2];
 char II[2];
 char IR[2];
 char typeRNA[10];
 char type[10];
 int iup;
 int idown;
 int exon;
 char one[40];
 char two[40];
 char three[40];
 char four[40];
 char ones[40];
 char twos[40];
 char threes[40];
 char fours[40];
 int gene_start;
 int gene_stop;
 char position[100];
 int in_start;
 int in_stop;
 int skip_start;
 int skip_stop;

 int right_intron_middle;
 int left_intron_middle;
 int exon_middle;

 int firstaregion_start;
 int firstaregion_stop;

 int secondaregion_start;
 int secondaregion_stop;

 int thirdaregion_start;
 int thirdaregion_stop;

 int fourtharegion_start;
 int fourtharegion_stop;

 double sumfirsta[30];
 double sumseconda[30];
 double sumthirda[30];
 double sumfourtha[30];

//////////////////

 int firstacheck[30];
 int secondacheck[30];
 int thirdacheck[30];
 int fourthacheck[30];

 double sumrowfirst;
 double sumrowsecond;
 double sumrowthird;
 double sumrowfourth;

 double newdirank;

 double max;

 int first_reg;
 int second_reg;
 int third_reg;
 int fourth_reg;

 double dIRank;
 
 double I1;
 double I2;
 double I3;
 double I4;
 double I5;
 double I6;
 double TC;
 double dI;

 double dIRankTDP;
 double dIRankTIA;
 char dIRankstr[40];
 double FDR;
 char FDRstr[40];
 char duplicate[20];
};
