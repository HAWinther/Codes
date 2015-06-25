
//=======================================================
// This file contains all the methods that relies on
// reading data from file and storing it....
//=======================================================

void init(int argv, char ** argc);
void read_halo_profile(std::string filename);
void read_halos(std::string filename);
int count_halos(std::string filename);

//=======================================================
// Read parameterfile with Lua
//=======================================================

void ReadFileWithLua::read(){
  int status;

  if(!isopen){
    std::cout << "LuaFile is not open so cannot read it!" << std::endl;
    exit(1);
  }

  // Read filenames etc.
  global.outfileprefix   = std::string(read_string2("outfileprefix",   &status, true));
  global.outfilesuffix   = std::string(read_string2("outfilesuffix",   &status, true));
  global.halofile        = std::string(read_string2("halofile",        &status, true));
  global.haloprofilefile = std::string(read_string2("haloprofilefile", &status, true));
  global.nbins           = read_int("nbins");
  global.boxsizeinkpc    = read_double("boxsizeinkpc");

  // Read profile mass-ranges
  global.nprofiles       = read_int("nprofiles");
  for(int i=0;i<global.nprofiles;i++){
    std::string tmp = "Mmin" + to_string(i+1);
    global.Mmin[i] = read_double(tmp.c_str());
  }
  for(int i=0;i<global.nprofiles;i++){
    std::string tmp = "Mmax" + to_string(i+1);
    global.Mmax[i] = read_double(tmp.c_str());
  }

  // Verbose
  std::cout << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << "Analyze halos running with...           "  << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << "Outfile prefix = " << global.outfileprefix << std::endl;
  std::cout << "Outfile suffix = " << global.outfilesuffix << std::endl;
  std::cout << "Halofile       = " << global.halofile      << std::endl;
  std::cout << "Haloprofile    = " << global.haloprofilefile << std::endl;
  std::cout << "Nbins          = " << global.nbins         << std::endl;
  std::cout << "Box (kpc/h)    = " << global.boxsizeinkpc  << std::endl;
  std::cout << "nprofiles      = " << global.nprofiles     << std::endl;
  for(int i=0;i<global.nprofiles;i++){
    std::cout << "Profile " << i+1 << " :" << global.Mmin[i] << " - " << global.Mmax[i] << " Msun/h" << std::endl;
  }
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
}

//=======================================================
// Get parameters from standard input
//=======================================================

void init(int argv, char ** argc){
  ReadFileWithLua paramfile;

  // Check if we have a filename
  if(argv==1){
    std::cout << "Run as ./analyzehalos param.lua" << std::endl;
    std::cout << "The format of the param.lua file is given in io.h (see ReadFileWithLua::read())" << std::endl;
    exit(1);
  } else {
    std::ifstream fp(argc[1]);
    if(!fp.good()){
      std::cout << "Parameterfile does not exist: " << argc[1] << std::endl;
      exit(1);
    }
    fp.close();
  }

  // Open and read parameters
  paramfile.open(argc[1]);
  paramfile.read();
  paramfile.close();
}

//=======================================================
// Read the halo file and store the halos in memory
//=======================================================

void read_halos(std::string filename){
  int npart, nhalos = global.nhalos, ncolhalofile = 43;
  realT x[NDIM], v[NDIM], rvir, mass, boxfac, tmp[ncolhalofile];
  Halo *halos = global.halos;
  std::string line;
  std::ifstream fp;

  std::cout << "=========================================" << std::endl;
  std::cout << "Read (AHF) Halo information" << std::endl;
  std::cout << "=========================================" << std::endl;

  fp.open(filename.c_str());
  if(!fp.good()){
    std::cout << "Halo file does not exist: " << filename << std::endl;  
    exit(1);
  }

  //==============================================================
  // ID(1)	hostHalo(2)	numSubStruct(3)	Mvir(4)	npart(5)	Xc(6)	
  // Yc(7)	Zc(8)	VXc(9)	VYc(10)	VZc(11)	Rvir(12) Rmax(13)	r2(14)
  // mbp_offset(15)	com_offset(16)	Vmax(17) v_esc(18)	sigV(19)
  // lambda(20)	lambdaE(21)	Lx(22)	Ly(23)	Lz(24) b(25)	c(26)	
  // Eax(27)	Eay(28)	Eaz(29)	Ebx(30)	Eby(31)	Ebz(32) Ecx(33)	
  // Ecy(34)	Ecz(35)	ovdens(36)	nbins(37)	fMhires(38) Ekin(39)	
  // Epot(40)	SurfP(41)	Phi0(42)	cNFW(43)	
  //==============================================================

  // Convert positions to code units [0,1]
  boxfac = 1.0/global.boxsizeinkpc;

  // Get first line
  getline(fp,line);

  // Rea AHF halo file
  for(int i=0;i<nhalos;i++){
    for(int j=0;j<ncolhalofile;j++)
      fp >> tmp[j];
    npart = int(tmp[4]);
    x[0]  = tmp[5]  * boxfac;
    x[1]  = tmp[6]  * boxfac;
    x[2]  = tmp[7]  * boxfac;
    rvir  = tmp[11] * boxfac;
    mass  = tmp[3];

    // If velocity is not given put v[i] to 0! Needs to be defined
    v[0] = tmp[8];
    v[1] = tmp[9];
    v[2] = tmp[10];

    // Create the halo and add it to the list
    halos[i] = Halo(x,v,rvir,mass,npart,i);
  }
  fp.close();

  // Verbose
  std::cout << "nhalos = " << nhalos << std::endl;
  std::cout << std::endl;
  halos[0].print_info();
  halos[global.nhalos-1].print_info();

  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
}

//=======================================================
// Count the number of lines in the halo file
// which gives us the number of halos + 1
//=======================================================

int count_halos(std::string filename){
  std::ifstream fp;
  int nhalos;
  fp.open(filename.c_str());
  if(!fp.good()){
    std::cout << "Error in count halos. Halofile file does not exist: " << filename << std::endl;  
    exit(1);
  }
  nhalos = std::count(std::istreambuf_iterator<char>(fp), 
      std::istreambuf_iterator<char>(), '\n') - 1;
  fp.close();
  return nhalos;
}

//=======================================================
// Read AHF profile output
// Assumes halos have been read and is in correct order
//=======================================================

void read_halo_profile(std::string filename){
  int ncolprofilefile = NCOLPROFILEFILE, nb, nbins[global.nhalos], i;
  realT tmp[ncolprofilefile], rold;
  std::ifstream fp;
  std::string line;

  //==========================================================
  // AHF profile-file column information:
  // r(1) npart(2)  M_in_r(3) ovdens(4) dens(5) vcirc(6)
  // vesc(7) sigv(8) Lx(9) Ly(10)  Lz(11)  b(12) c(13) Eax(14)
  // Eay(15) Eaz(16) Ebx(17) Eby(18) Ebz(19) Ecx(20) Ecy(21)
  // Ecz(22) Ekin(23)  Epot(24)
  //==========================================================

  std::cout << "=========================================" << std::endl;
  std::cout << "Read Halo Profiles" << std::endl;

  fp.open(filename.c_str());
  if(!fp.good()){
    std::cout << "Profile file does not exist!" << filename << std::endl;  
    exit(1);
  }

  // Read first line
  getline(fp,line);

  // Read through file and get number of bins for each halo
  int totlines = 0;
  rold = 1e30;
  nb = i = 0;
  while(1){
    fp >> tmp[0];
    if(fp.eof()){
      nbins[i++] = nb;
      break;
    }

    // Check if new profile has started...
    if(tmp[0] < 0.0 && rold > 0.0 && nb>1){
      nbins[i++] = nb;
      //std::cout << "Profile " << i << " nbins = " << nb << std::endl;
      nb = 0;
    }
    rold = tmp[0];

    // Sanity check
    if(i>global.nhalos){
      std::cout << "Error: Too many profiles found in the halofile ";
      std::cout << "i= "<< i << " nhalos= " << global.nhalos << std::endl;
      exit(1);
    }
    if(nb > 50){
      std::cout << "Error: Too many bins in profile? nb>50 i = " << i << " totlinesread = " << totlines << std::endl;
      exit(1);
    }

    // Read rest of column
    getline(fp,line);
    ++nb;
    ++totlines;
  }

  std::cout << "Done read ... rewind and add!" << std::endl;

  // Rewind the file and read first line
  fp.clear();
  fp.seekg(0);
  getline(fp,line);

  // Read through it again and store data
  for(i=0;i<global.nhalos;i++){
    // Allocate memory for the profile
    global.halos[i].allocate_memory_profile(nbins[i], ncolprofilefile-1);

    for(int j=0;j<nbins[i];j++){
      // Read data
      for(int k=0;k<ncolprofilefile;k++) fp >> tmp[k];

      // Store the data
      global.halos[i].store_profile_data(tmp,j);
    }
  }

  std::cout << "Spline all profiles" << std::endl;
  for(i=0;i<global.nhalos;i++){
    global.halos[i].spline_halo_profiles();
  }

  std::cout << "Done..." << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
}

