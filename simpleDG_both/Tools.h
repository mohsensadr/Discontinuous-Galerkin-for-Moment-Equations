#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <iostream>
#include <algorithm>

using namespace std;

// Store the formatted string of time in the output
void GetTimeStamp( char *output );

// simple time measuring tools
void startTimer( string s );
void stopTimer();
double GlobalTime();

// command line parameters
int exists_argument(int argc,const char **argv,const char *keystr);
float get_float_argument(int argc,const char **argv,const char *keystr);
int get_string_argument(int argc,const char **argv,const char *keystr,char *result);

// sorting
struct Ordering
{
  Ordering( VectorXd &Array );
  bool operator() (int i,int j);    
  VectorXi index;  
  VectorXd array;
};


