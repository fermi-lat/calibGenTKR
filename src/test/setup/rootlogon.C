{
    char* macros[]={
        "RootTreeAnalysis.cxx",
        "NtupleAnalysis.cxx",
        "chainTrees.cxx",
        "chainAll.cxx",
        "copyTree.cxx",
        "createEventList.cxx",
        "pruneTree.cxx"};
    char* libs[]={
        "mcRootData",
        "digiRootData",
        "reconRootData"};

    printf("\nLoading GLAST macros...\n");
    for(int i=0; i< sizeof(macros)/sizeof(char*); i++){
        printf("\t%s\n", macros[i]);
        gROOT->LoadMacro(macros[i]);
    }
    printf("\nLoading GLAST libraries...\n");
    // Load the appropriate version - depending on what OS 
    // we are running
    if (!strcmp(gSystem->GetName(), "WinNT")) {
        // nt version can use CMT for env vars.
        char buf[128]; 
        for(int j=0; j< sizeof(libs)/sizeof(char*); j++){
           char * env = getenv( strcat( strcpy(buf,libs[j]), "Shr") );
           // If we are not using CMT - we will use the ROOTANALYSIS env var to determine
           // where the GLAST ROOT libraries are located.
           char * raPath = getenv("ROOTANALYSIS");
           char * path = env? strcat( strcpy(buf,env), ".dll")
                            : strcat( strcat( strcat( strcpy(buf, raPath), "/lib/") ,libs[j]),".dll");
           printf("\t%s\n",path);
           gSystem->Load( path ); 
        }
       
    } else {  // UNIX
        gSystem->Load("libPhysics.so");
        gSystem->Load("libmcRootData.so");
        gSystem->Load("libdigiRootData.so");
        gSystem->Load("libreconRootData.so");
    }
    
}

