//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


/*
     Comparion functions for floating point types (FPT)
     
     Oliver Fochler

     last modified: 24.10.2006
*/


#ifndef FPT_COMPARE
#define FPT_COMPARE


#define FPT_COMP_PRECISION 1.0e-3

bool FPT_COMP_E(const double a, const double b);  // check for a == b
bool FPT_COMP_E(const float a, const float b);


bool FPT_COMP_L(const double a, const double b);  // check for a < b
bool FPT_COMP_L(const float a, const float b);

bool FPT_COMP_G(const double a, const double b);  // check for a > b
bool FPT_COMP_G(const float a, const float b);

bool FPT_COMP_LE(const double a, const double b); // check for a <= b
bool FPT_COMP_LE(const float a, const float b);

bool FPT_COMP_GE(const double a, const double b); // check for a >= b
bool FPT_COMP_GE(const float a, const float b);



#endif  //#ifndef FPT_COMPARE
