
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
 
	 ret = tetramer(opt);
	}catch (exception& e){
		cerr << "exception caught: " << e.what() << endl;
	}
return 0;
}
