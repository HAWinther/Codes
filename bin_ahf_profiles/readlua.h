#ifndef READLUAHEADER_INC
#define READLUAHEADER_INC

extern "C"{
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

//=======================================================
// Read a file with Lua class
//=======================================================

class ReadFileWithLua{
  private:
    lua_State* L;
    std::string filename;
    bool isopen;
  public:

    ReadFileWithLua(){
      isopen=false;
    }

    ReadFileWithLua(std::string filename){
      this->filename = filename;
      open(filename);
      isopen = true;
    }

    ~ReadFileWithLua(){
      if(isopen) lua_close(L);
    }

    void open(std::string filename){
      this->filename = filename;
      L = luaL_newstate();
      luaL_openlibs(L);
      if(luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0)) {
        std::cout << "Error in reading parameterfile: " << lua_tostring(L, -1) << std::endl;
        exit(1);
      }
      isopen=true;
    }

    void read();

    void close(){
      if(isopen){
        lua_close(L);
        isopen = false;
      } else {
        std::cout << "LuaFile is already closed..." << std::endl;
      }
    }

    int read_int(const char name[]){
      lua_getglobal(L, name);
      if(!lua_isnumber(L, -1)) {
        std::cout << "Parameter " << name << " not found in inputfile" << std::endl;
        exit(1);
      }
      int n = lua_tointeger(L, -1);
      lua_pop(L, 1);
      return n;
    }

    double read_double(const char name[]){
      lua_getglobal(L, name);
      if(!lua_isnumber(L, -1)) {
        std::cout << "Parameter " << name << " not found in inputfile" << std::endl;
        exit(1);
      }
      double val = lua_tonumber(L, -1);
      lua_pop(L, 1);
      return val;
    }

    char* read_string2(const char name[], int* len, bool required){
      lua_getglobal(L, name);
      if(!lua_isstring(L, -1)) {
        if(required){
          std::cout << "Error: Parameter " << name << " not found in the parameter file" << std::endl;
          exit(1);
        } else {
          *len = 0;
          return 0;
        }
      }
      char const * const str = lua_tostring(L, -1);
      const int n = strlen(str) + 1;
      char* const val =  new char[n];
      strncpy(val, str, n);
      lua_pop(L, 1);
      *len = n;
      return val;
    }

    bool read_bool(const char name[]){
      lua_getglobal(L, name);
      if(!lua_isboolean(L, -1)) {
        std::cout << "Error: Parameter %s not found in the parameter file" << name << std::endl;
      }
      int n = lua_toboolean(L, -1);
      lua_pop(L, 1);
      return bool(n);
    }

    double* read_array_double(const char name[], int *len){
      int i;
      lua_getglobal(L, name);
      if(!lua_istable(L, -1)) {
        std::cout << "Error: Parameter %s not found in the parameter file" << name << std::endl;
      }
      const int n = luaL_len(L, -1);
      double* const array = new double[n];
      for(i=1; i<=n; ++i) {
        lua_pushinteger(L, i);
        lua_gettable(L, -2);
        double x = lua_tonumber(L, -1);
        lua_pop(L,1);
        array[i-1] = x;
      }
      lua_pop(L, 1);
      *len = n;
      return array;
    }
};

#endif
