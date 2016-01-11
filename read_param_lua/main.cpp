#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "readlua.h"

//=======================================================
// Code to read parameterfiles in LUA format
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//=======================================================

//=======================================================
// Implementation of the read function
//=======================================================

void ReadFileWithLua::read(){
  std::string inpstring;
  std::vector<double> inparray;
  int inpint;
  double inpdouble;
  bool status;

  // Check if file is open
  if(!isopen){
    std::cout << "LuaFile is not open so cannot read it!" << std::endl;
    return;
  }

  // Read different parameters from lua file. False = we don't abort it param is not found

  // Look for a string named 'inpstring'
  inpstring  = read_string      ("inpstring", &status, false);
  if(status) std::cout << "String 'inpstring' found : " << inpstring << std::endl;
  std::cout << std::endl; 

  // Look for an integer named 'inpint'
  inpint     = read_int         ("inpint",    &status, false);
  if(status) std::cout << "Integer 'inpint'   found  : " << inpint << std::endl;
  std::cout << std::endl; 

  // Look for an double named 'inpdouble'
  inpdouble  = read_double      ("inpdouble", &status, false);
  if(status) std::cout << "Double 'inpdouble' found  : " << inpdouble << std::endl;
  std::cout << std::endl; 

  // Look for an double array named 'inparray'
  inparray   = read_array<double>("inparray",  &status, false);
  if(status){
    std::cout << "Array 'inparray'   found  : ";
    for(int i=0;i<inparray.size();i++)
      std::cout << (i==0 ? "" : ", ") << inparray[i];
    std::cout << std::endl; 
  }
}

int main(int argv, char **argc){
  std::string filename;
  if(argv == 1) {
    std::cout << "Parameterfile not provided.\nRun as ./readlua input.lua" << std::endl;
    exit(1);
  }
  filename = argc[1];

  // Open and read parameters
  ReadFileWithLua paramfile;
  paramfile.open_read_close(filename);

  return 0;
}

