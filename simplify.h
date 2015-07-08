#ifndef simplify_h
#define simplify_h 1



#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>



typedef  double G4double;  // for consistency
typedef int G4int;

using namespace std;

typedef struct mypoint {
     double e;
     double xs; } Point;

int  simplify_function(G4double ,
		       vector <Point> & ,
		       vector<Point>  &);

void RemoveBias( vector <Point> &,
		 vector <Point> &,
		 vector <Point> &);


#endif
