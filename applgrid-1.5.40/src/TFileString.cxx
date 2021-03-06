/**
 **     @file    TFileString.cxx
 **
 **     @brief   root TObject std::string std::vector class for writing std::string std::vectors
 **              to root files               
 **
 **     @author  mark sutton
 **     @date    Sat Mar 15 19:50:15 GMT 2008 
 **
 **     @copyright (C) 2002-2019 mark sutton (sutt @ cern.ch) 
 **
 **     $Id: TFileString.cxx, v0.0   Sat Mar 15 19:50:15 GMT 2008 sutt $
 **
 **/



#include "TFileString.h"

ClassImp(TFileString)


std::ostream& operator<<(std::ostream& s, const TFileString& fs) { 
  s << fs.name() << ":";
  for ( unsigned i=0 ; i<fs.size() ; i++ ) s << "\t" << fs[i];
  return s;
}



