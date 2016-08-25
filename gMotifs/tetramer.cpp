
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "inc.h"

int main(int argc, char *argv[]){
 try{
  cout << endl;
  /*-----------------------------------------------
                  CONFIGURATION
   ------------------------------------------------*/
	 options opt(argc,argv);
	 opt.print(cout);

	 int ret = 0;
 
	 if( argc<2 ){
		 cout << "Missing splicing change filename" << endl;
		 return -1;
	 } 

	 string splicing_change(argv[1]),
  	 percent = string("percent_0.5/");
	 ret = tetramer(opt,splicing_change,percent);
	}catch (exception& e){
		cerr << "exception caught: " << e.what() << endl;
	}
return 0;
}
