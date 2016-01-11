#ifndef READFILEWITHLUA_INC
#define READFILEWITHLUA_INC

#include <stdlib.h>
#include <string>
#include <vector>

extern "C"{
  #include <lua.h>
  #include <lauxlib.h>
  #include <lualib.h>
}

//=======================================================
//
// Methods to read parameterfiles Lua style
//
// Allows for paramfiles like in thie code, e.g.:
//   vara = 2.0
//   varb = {1.0, 2.0, 5.0}
//   varc = "Hello"
//
// For arrays only integer and floating point types are
// implemented (so no string and bool arrays yet).
//
//=======================================================

class ReadFileWithLua{
  private:
    lua_State* L;
    std::string filename;
    bool isopen;
  public:

    ReadFileWithLua() : isopen(false) {}

    ReadFileWithLua(std::string filename){
      open(filename);
    }

    ~ReadFileWithLua(){
      if(isopen) lua_close(L);
    }

    //====================================================
    // Open file
    //====================================================

    void open(std::string filename){
      this->filename = filename;
      L = luaL_newstate();
      luaL_openlibs(L);
      if(luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0)) {
        std::cout << "Error in reading parameterfile: " << lua_tostring(L, -1) << std::endl;
        isopen = false;
      }
      isopen = true;
    }

    void open_read_close(std::string filename){
      open(filename);
      read();
      close();
    }

    //====================================================
    // Close the file
    //====================================================

    void close(){
      if(isopen){
        lua_close(L);
      } else {
        std::cout << "LuaFile " << filename << " is already closed!" << std::endl;
      }
      isopen = false;
      L = NULL;
    }

    //====================================================
    // Read function for user to implement
    //====================================================

    void read();

    //====================================================
    // Read functions. Should template these functions
    // but it's a mess so one function per type...
    //====================================================

    //====================================================
    // Read a single integer
    //====================================================

    int read_int(std::string name, bool *found, bool required = true){
      int val = 0;

      // Check if value is in file
      lua_getglobal(L, name.c_str());
      *found = lua_isnumber(L, -1);
      if(!*found) {
        if(required){
          std::cout << "Parameter " << name << " not found in inputfile" << std::endl;
          exit(1);
        }
        return val;
      }

      // Fetch value
      val = lua_tointeger(L, -1);
      lua_pop(L, 1);
      return val;
    }

    //====================================================
    // Read a double
    //====================================================

    double read_double(std::string name, bool *found, bool required = true){
      double val = 0.0;

      // Check if value is in file
      lua_getglobal(L, name.c_str());
      *found = lua_isnumber(L, -1);
      if(!*found) {
        if(required){
          std::cout << "Parameter " << name << " not found in inputfile" << std::endl;
          exit(1);
        }
        return val;
      }

      // Fetch value
      val = lua_tonumber(L, -1);
      lua_pop(L, 1);
      return val;
    }

    //====================================================
    // Read a string
    //====================================================

    std::string read_string(std::string name, bool *found, bool required = true){
      std::string val = "";

      // Check if parameter is in file
      lua_getglobal(L, name.c_str());
      *found = lua_isstring(L, -1);
      if(!*found) {
        if(required){
          std::cout << "Error: Parameter " << name << " not found in the parameter file" << std::endl;
          exit(1);
        }
        return val;
      }
     
      // Fetch value
      val = lua_tostring(L, -1);
      lua_pop(L, 1);
      return val;
    }

    //====================================================
    // Read a single bool
    //====================================================

    bool read_bool(std::string name, bool *found, bool required = true){
      bool val = false;

      // Check if parameter is in file
      lua_getglobal(L, name.c_str());
      *found = lua_isboolean(L, -1);
      if(!*found) {
        if(required){
          std::cout << "Error: Parameter " << name << " not found in the parameter file" << std::endl;
          exit(1);
        }
        return val;
      }

      // Fetch value
      val = bool(lua_toboolean(L, -1));
      lua_pop(L, 1);
      return val;
    }

    //====================================================
    // Read an array (only works for int/double types as
    // written)
    //====================================================

    template <class T>
    std::vector<T> read_array(std::string name, bool *found, bool required = true){
      int n;
      std::vector<T> val;

      // Check if parameter is in file
      lua_getglobal(L, name.c_str());
      *found = lua_istable(L, -1);
      if(!*found) {
        if(required){
          std::cout << "Error: Parameter " << name << " not found in the parameter file" << std::endl;
          exit(1);
        }
        return val;
      }

      // Fetch value(s)
      n = luaL_len(L, -1);
      for(int i=1;i<=n;i++){
        lua_pushinteger(L, i);
        lua_gettable(L, -2);
        val.push_back(T(lua_tonumber(L, -1)));
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
      return val;
    }
};
#endif
