void analysis();

//=======================================================
// Driver routine for the whole analysis
//=======================================================

void analysis(){
  int nbins = global.nbins;

  //==============================================
  // Halos
  //==============================================

  // Create halos
  global.nhalos = count_halos(global.halofile);
  Halo *halos   = new Halo[global.nhalos];
  global.halos  = halos;

  // Read halos from file
  read_halos(global.halofile);

  // Read halo profiles from file
  read_halo_profile(global.haloprofilefile);

  // Output binned profiles from AHF data
  for(int i=0;i<global.nprofiles;i++)
    bin_ahf_profile(global.Mmin[i], global.Mmax[i], nbins, global.outfileprefix + to_string(i+1) + global.outfilesuffix);

  // Clean up
  delete[] halos;
}

