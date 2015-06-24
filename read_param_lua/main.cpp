#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "readlua.h"

/////////////////////////////////////////
// Implementation of the read function
/////////////////////////////////////////

void ReadFileWithLua::read(){
  int status;
  std::string mystring;
  int myint;
  double mydouble;

  if(!isopen){
    std::cout << "LuaFile is not open so cannot read it!" << std::endl;
    exit(1);
  }

	// Read from lua file
	mystring  = std::string(read_string2("mystring", &status, true));
	myint     = read_int("myint");
	mydouble  = read_double("mydouble");

	// Output the data read from the lua file
	std::cout << "String : " << mystring << std::endl;
	std::cout << "Int    : " << myint    << std::endl;
	std::cout << "Double : " << mydouble << std::endl;
}

/////////////////////////////////////////
// Open a lua paramfile and read it
/////////////////////////////////////////

int main(int argv, char **argc){
	ReadFileWithLua paramfile;

  if(argv == 1) exit(1);

	// Open and read parameters
	paramfile.open(argc[1]);
	paramfile.read();
	paramfile.close();

	return 0;
}

