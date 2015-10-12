{
  // Set up FW Lite for automatic loading of CMS libraries
  // and data formats.   As you may have other user-defined setup
  // in your rootlogon.C, the CMS setup is executed only if the CMS
  // environment is set up.
  //
  //TString cmsswbase = getenv("CMSSW_BASE");
  //if (cmsswbase.Length() > 0) {
  //  //
  //  // The CMSSW environment is defined (this is true even for FW Lite)
  //  // so set up the rest.
  //  //
  //  cout << "Loading FW Lite setup." << endl;
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  
  gSystem->Load("libPhysicsToolsTagAndProbe.so");
  //  //gSystem->Load("libDataFormatsPatCandidates.so");
  //}
  //gSystem->Load("utilities/RooCBExGaussShape_cxx.so");
}
