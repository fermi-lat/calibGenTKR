/** @file badStripsCalib.cxx
@brief Algorithm to calibrate bad strips
*/

#include "badStripsCalib.h"
#include <string>
#include <algorithm>
#include "facilities/commonUtilities.h"
UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum,  reconRunNum,  mcRunNum;

static Int_t evtCount = 0;

// some test flags
bool overlayTowers   = false;
bool debug      = false;
bool debugf     = false;
bool infof      = false;
bool readOnly   = false; // for timing the overhead... just read the events, no computation

// choose fit method
bool simpleCuts   = true;   // dead = no hits; hot = occupancy > maxOccupancy (above average for plane)
bool testAgainstFit = false;  // code to handle the edges, where there are few counts

// used for average or local-average cuts
const float nSig = 6.f;
const float hotFactor = 10.f;
//const float maxOccupancy = 0.01f; -> now a member, setable in options file

// checkCount: if true, require an average of at least minHistCount counts before testing for dead
// not completely thought out, but shouldn't make much difference in a flight tower
const bool  checkCount   = false;
const int   minHistCount = 16;

const int laddersPerPlane = 4;
const int stripsPerChip   = 64;
const int stripsPerBlock  = stripsPerChip/2;
const int chipsPerLadder  = 6;
const int stripsPerLadder = stripsPerChip*chipsPerLadder;
const int stripsPerPlane  = stripsPerLadder*laddersPerPlane;
const int numChips = chipsPerLadder*laddersPerPlane;
const int numBlocks = numChips*2;

//const int numHists = numPlanes*numTowers;
int numLadders[36]; // historical from balloon flight (not likely to ever need this!)

int i,j;

// Some special histograms
TH1F *outsideTowers, *outsidePlanes;

void BadStripsCalib::DigiHistDefine() {
    // Purpose and Method:  Digitization histogram definitions

    // Must cd into the histFile to be sure these new histograms are created in
    // this histogram ROOT file
    histFile->cd();

    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111);

    char buf[10];
    char buf1[40];

    for(i=0;i<m_nPlanes;++i) {numLadders[i] = laddersPerPlane;}

    int nTowers = m_towerNums.size();
    int nHists = nTowers*m_nPlanes;

    m_tkrHists.assign(nHists, static_cast<TH1F*>(NULL));

    for (j = 0; j<nTowers; ++j) {
        int tower = m_towerNums[j];
        for (i = 0; i<m_nPlanes; ++i) {
            sprintf(buf,"TKR %02d-%02d",tower,i);
            sprintf(buf1,"Strips hit in tower %d, plane %d",tower, i);
            int nStrips = numLadders[i]*stripsPerLadder;
            m_tkrHists[m_nPlanes*j +i] = new TH1F(buf, buf1, nStrips, -0.5, nStrips-0.5);
        }
    }

    outsideTowers = new TH1F("outsideTowers", "hits in Towers outside current config", 18, -1.5, 16.5);
    outsidePlanes = new TH1F("outsidePlanes", "hits in planes outside current config", 38, -1.5, 36.5);
    return;       
}

void BadStripsCalib::DigiTkr() {
    // Purpose and Method: Process on TKR digi event

    const TObjArray* tkrDigiCol = evt->getTkrDigiCol();
    if (!tkrDigiCol) return;

    int nTkrDigi = tkrDigiCol->GetEntries();
    if (debug) std::cout << "Number of Digis: " << nTkrDigi << std::endl;

    TIter tkrDigiIter(tkrDigiCol);
    TkrDigi *t = 0;

    //int numTowers = m_towerNums.size();

        // yes, it's "=", not "=="... we're assigning the pointer to the next object to t...
        while (t = (TkrDigi*)tkrDigiIter.Next()) {

        int tower = (t->getTower()).id();
        if (overlayTowers) tower = 0;
        if( tower>=16 || tower < 0) {
            int truncTower = std::max(-1, std::min(16, tower));
            outsideTowers->Fill(truncTower);
            continue;
        }
        int layer = t->getBilayer();
        GlastAxis::axis view = t->getView();
        int iview = (int) view;
        bool isY = (view==1);
        int plane = 2*layer + (1-isY)*(1-layer%2) + isY*layer%2;
        if (plane>=m_nPlanes || plane < 0) {
            int truncPlane = std::max(-1, std::min(36, plane));
            outsidePlanes->Fill(truncPlane);
            continue;
        } 
        if (debug) std::cout << "    tower,layer,view " << tower << " " 
            << layer << " " << iview <<std::endl;

        unsigned int numStrips = t->getNumHits();
        unsigned int iHit;
        for (iHit = 0; iHit < numStrips; ++iHit) {
            int strip = t->getHit(iHit);
            if (debug) std::cout << strip << " " ;
            int histId = m_histId[tower];
            if (histId>-1) {
                m_tkrHists[m_histId[tower]*m_nPlanes + plane]->Fill(strip);
            } else {
                outsideTowers->Fill(tower);
            }
        }
    }
    return;
}

void BadStripsCalib::HistDefine() {
    // Purpose and Method:  Setup Histograms

    histFile = new TFile(m_histFileName,"RECREATE");

    //McHistDefine();
    DigiHistDefine();
    //ReconHistDefine();
}

void BadStripsCalib::Go(Int_t numEvents)
{    
    // Purpose and Method:  Event Loop
    //   All analysis goes here

    //  To read only selected branches - saves processing time
    //  Comment out any branches you are not interested in.

    if (digiTree) {
        digiTree->SetBranchStatus("*",0);  // disable all branches
        // activate desired branches
        //digiTree->SetBranchStatus("m_cal*",1);  
        digiTree->SetBranchStatus("m_tkr*",1);  
        //digiTree->SetBranchStatus("m_acd*",1);
        digiTree->SetBranchStatus("m_eventId", 1); 
        digiTree->SetBranchStatus("m_runId", 1);
    }

    if (m_digiChain) {
        m_digiChain->SetBranchStatus("*",0);  // disable all branches
        // activate desired branches
        //m_digiChain->SetBranchStatus("m_cal*",1);  
        m_digiChain->SetBranchStatus("m_tkr*",1);  
        //m_digiChain->SetBranchStatus("m_acd*",1);
        m_digiChain->SetBranchStatus("m_eventId", 1); 
        m_digiChain->SetBranchStatus("m_runId", 1);
    }

    if (reconTree) {
        reconTree->SetBranchStatus("*",0);
    }
    if (m_recChain) {
        m_recChain->SetBranchStatus("*",0);
    }

    // determine how many events to process
    Int_t nentries = GetEntries();
    std::cout << "\nNum Events in File is: " << nentries << std::endl;
    Int_t curI;
    Int_t nMax = TMath::Min(numEvents+m_StartEvent,nentries);
    Int_t nEventsToRead = TMath::Min(numEvents+m_StartEvent,nentries-m_StartEvent);
    std::cout << nEventsToRead << " events will be read, starting at event " 
              << m_StartEvent << std::endl;

    if (m_StartEvent == nentries) {
        std::cout << " all events in file read" << std::endl;
        return;
    }
    if (nentries <= 0) return;

    if(readOnly) std::cout << "No event processing will be done" << std::endl;

    // Keep track of how many bytes we have read in from the data files
    Int_t nbytes = 0, nb = 0;

    TStopwatch timer, timer1;
    timer.Start();
    timer1.Start();


    // BEGINNING OF EVENT LOOP
    int nTimeCheck;
    if      (nEventsToRead<50000)  {nTimeCheck = 5000;}
    else if (nEventsToRead<100000) {nTimeCheck = 10000;}
    else if (nEventsToRead<200000) {nTimeCheck = 20000;}
    else                           {nTimeCheck = 50000;}
    for (Int_t ievent=m_StartEvent; ievent<nMax; ievent++, curI=ievent) {
        int nProcessed = ievent-m_StartEvent;
        if (nProcessed%nTimeCheck==0 && nProcessed>0) {
            std::cout << ievent-m_StartEvent << " events processed, " 
                << " Real time " << timer.RealTime()
                << " Cpu time " << timer.CpuTime() << std::endl;
            timer.Start();
        }

        evtCount++;

        //if (mc) mc->Clear();
        if (rec) rec->Clear();

        if (evt) evt->Clear();

        digiEventId = 0; reconEventId = 0; mcEventId = 0;
        digiRunNum = 0; reconRunNum = 0; mcRunNum = 0;


        nb = GetEvent(ievent);
        nbytes += nb;

        if (readOnly) continue;

        // Digi ONLY analysis
        if (evt) {
            digiEventId = evt->getEventId(); 
            digiRunNum = evt->getRunId();
            DigiTkr();
            //DigiCal();
            //DigiAcd();
        }    
        if (rec) {
            reconEventId = rec->getEventId();
            reconRunNum  = rec->getRunId();
        }
    }  // end analysis code in event loop

    m_StartEvent = curI;

    std::cout << "total time: Real " << timer1.RealTime()
        << " Cpu " << timer1.CpuTime() << std::endl;
}

/** @todo
  figure out better local average, for now we use simple average
*/

void BadStripsCalib::Finish()
{
    std::cout << std::endl << evtCount << " events processed" << std::endl << std::endl;
    
    if(readOnly) return;

    // Get the time stamp of first and last event

    double startTime, stopTime;
    std::string time1, time2;
    GetEvent(0);
    if (evt) {
        digiEventId = evt->getEventId(); 
        digiRunNum = evt->getRunId();
        startTime = evt->getTimeStamp();
        // for now
        time1 = "2003-10-1-00:00";       
        evt->Print();
        std::cout << "First event: Run, event, time " << digiRunNum << " " 
            << digiEventId << " " << startTime << " " << time1.c_str() <<  std::endl;
    }

    int nEvts = GetEntries();
    if (infof) std::cout << "nEvts" << nEvts << std::endl;
    GetEvent(nEvts-1);
    if (evt) {
        digiEventId = evt->getEventId(); 
        digiRunNum = evt->getRunId();
        stopTime = evt->getTimeStamp();
        // for now
        time2 = "2003-12-1-00:00";
        evt->Print();
        std::cout << "First event: Run, event, time " << digiRunNum << " " 
            << digiEventId << " " << stopTime << " " << time2.c_str() << std::endl;
    }

    int ihist;
    int numHists = m_tkrHists.size();
    std::vector<float> average(numHists);
    std::vector<bool>  towerCount(numHists);
    std::vector< std::vector<int> > chipCount(numHists, vector<int>(2*numChips));

    histFile->cd();

    int nTowers = m_towerNums.size();
    for (i=0;i<nTowers; ++i) { 
        //TKRCHIP[i]->Reset();
        towerCount[i] = false;
    }

    if (infof) std::cout << std::endl << "Analyzing strip histos" << std::endl;
    for (ihist=0; ihist<numHists; ++ihist){
        average[ihist] = 0;
        for (i=0;i<numChips*2;++i) { chipCount[ihist][i] = 0; }

        int isum = 0;
        double binsum = 0;
        int plane = ihist%m_nPlanes;
        int iTower = ihist/m_nPlanes;
        int tower = m_towerNums[iTower];
        if (debugf)
            std::cout << "tower, plane " << tower << " " << plane << std::endl;
        int nbin = numLadders[plane]*stripsPerLadder;
        int bin = nbin;

        // get rough average
        int entries = (int) m_tkrHists[ihist]->GetEntries();
        
        // refine the average a bit, by dropping the very low and very high bins
        if(entries>0) {
            towerCount[tower] = true;
            float xave = entries;
            xave /= nbin;
            if (debugf) 
                std::cout << " entries, ave " << entries << " " << xave << std::endl;

            for (bin=1; bin<nbin+1; ++bin) {
                float xbin = m_tkrHists[ihist]->GetBinContent(bin);
                //sumhist->Fill(xbin/xave);

                int block = (bin-1)/stripsPerBlock;
                chipCount[ihist][block] += (int) xbin;
                if ((xbin>0.1*xave && xbin<3*xave) /*|| (xbin>0 && xbin<10)*/) {
                    isum++; binsum += xbin;
                }
                if (isum==0) continue;
            }
            average[ihist] = binsum/isum;
        }

        if (debugf) std::cout << binsum << " counts for "<<isum
                         <<" bins, average = " << average[ihist] << std::endl;
        if (debugf) {
            std::cout << "chipCount, hist" << ihist << ": ";
            for (i=0;i<numChips*2;++i) {
                std::cout << chipCount[ihist][i] << " " ;
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl << std::endl;

    if (testAgainstFit) {
        for (ihist=0; ihist<numHists; ++ihist){
            std::cout << "Fitting m_tkrHists " << ihist << std::endl;
            if (average[ihist]>0) m_tkrHists[ihist]->Fit("pol8");    
        }
    }

    double xval[1];
    double yval;

    int badCount;
    bool testCount;

    //const int ntest = 2;
    float cutval[2] = {0, 0};

    // Dead strips 
    // call it dead if it is more than n sigma below the mean (or zero if that is greater)
    //    nSig is currently set to 6 which gives a ~10-7 chance of false dead strip;
    //    or ~0.1 per calibration

    // Hot strips 
    // for purposes of the offline calibration, call it hot if it fires more than 
    //    maxOccupancy of the time (over the average) or nSig sigma greater than the mean,
    //    if this is greater...

    int itest;
    std::string stype[2] = {"dead", "hot"};
    std::string calType[2] = {"TkrBadStrips", "TkrBadStrips"};
    std::string outFile[2] = {"deadStrips.xml", "hotStrips.xml"};
    std::ofstream listout;
    std::string indent1 = "   ", indent2 = "      ", indent3 = "         ",
        indent4 = "            ";
    std::string howBad = " nOnbdTrig=\"false\" nOnbdCalib=\"false\" nOnbdData=\"true\" ";

    for (itest=0; itest<2;++itest) {
        int inTower = -1;
        bool firstTower = true;
        bool flagIt;
        
        std::string outputFile;
        outputFile = m_xmlPath+m_prefix+"_"+outFile[itest];
        std::cout << "Output file " << outputFile.c_str() << std::endl;
        listout.open(outputFile.c_str(), std::ios::out);

        std::string buffer;
        std::ifstream infile;

        //std::string path = ::getenv("CALIBGENTKRROOT");
	std::string dtdPath = facilities::commonUtilities::joinPath(facilities::commonUtilities::getXmlPath("CalibGenTKR"), "badStrips.dtd");
	//std::string dtdPath = path+"/xml/badStrips.dtd";

        infile.open(dtdPath.c_str(), std::ios::in | std::ios::binary);
        if (!infile.is_open()) { std::cout << "dtd file not found!" << std::endl;
        return;
        }
        while (!infile.eof()) {
            getline(infile,buffer);
            if (debugf) std::cout << "*"<<  buffer << "*" <<std::endl;
            listout << buffer;
        }
        infile.close();

        if (debugf) std::cout << "test " << itest << std::endl;

        // write out the generic stuff: need to get current date
        // and maybe some of the other hard-wired stuff

        // here we need to get the current date
        std::string timeStamp = "2003-11-7-18:00";

        listout << std::endl <<
            "<badStrips badType= \"" << stype[itest].c_str() << "\">" << std::endl;

        listout << indent1.c_str()
            << "<generic instrument=\"EM\" timestamp=\"" << timeStamp.c_str() << "\"" 
            << " calType=\"" << calType[itest].c_str() << "\" fmtVersion=\"v2r0\" >" 
            << std::endl;


        listout << indent2.c_str() << "<inputSample startTime=\"" << time1.c_str()   
            << "\" stopTime=\"" << time2.c_str() 
            << "\" triggers=\"physics\" mode=\"normal\" source=\"XXX\" >" 
            << std::endl;

        listout << indent3.c_str() << "Output from BadStripsCalib"
            << std::endl;
        listout << indent2.c_str() << "</inputSample>" << std::endl;
        listout << indent1.c_str() << "</generic>" << std::endl;

        for (ihist=0; ihist<numHists; ++ihist){
            badCount = 0; testCount = false;
            int totalBadCount = 0;
            bool newHist = true;
            int iTower = ihist/m_nPlanes;
            int tower = m_towerNums[iTower];
            int towerCol = tower/4;
            int towerRow = tower - 4*towerCol;

            // tower is all bad
            if (!towerCount[tower] && itest==0) {
                if (inTower!=tower) {
                    if(!firstTower) listout << indent1.c_str() << "</tower>" << std::endl;
                    listout << indent1.c_str() << "<tower row=\"" << towerRow 
                        << "\" col=\"" << towerCol << "\""
                        << " allBad=\"true\"" << howBad << " >" << std::endl;
                    inTower = tower;
                    firstTower = false;
                }
                continue;
            }

            int plane = ihist%m_nPlanes;
            int nbin = numLadders[plane]*stripsPerLadder;

            int tray = (plane+1)/2;
            int botTop = (plane+1)%2;

            if (simpleCuts) {
                cutval[0] = 0;
                // call a strip hot if it fires more than maxOccupancy of the time (over the average)
                cutval[1] = average[ihist] + evtCount*m_maxOccupancy; 
            } else if (!testAgainstFit) {
                cutval[0] = std::max((float)0.0,(float)(average[ihist] - nSig*sqrt(average[ihist])));
                cutval[1] = std::max((float)(hotFactor*average[ihist]), 
                    (float)(average[ihist] + nSig*sqrt(average[ihist])));
            }

            float ave1 = average[ihist];
            if(infof) {
                std::cout << "Average count: " << ave1 << " nSig: " << nSig 
                    << " hotFactor: " << hotFactor ;
                std::cout << " cutVals: " << cutval[0] << " " << cutval[1] << std::endl;
            }

            int bin;
            TF1* fitfunc;
            if(testAgainstFit) fitfunc = m_tkrHists[ihist]->GetFunction("pol8");

            //float xval0 = 0;
            if (itest==1 || (average[ihist]>minHistCount || average[ihist]==0) ) {
                for (bin=1; bin<=nbin+1; ++bin) {
                    float xbin = m_tkrHists[ihist]->GetBinContent(bin);
                    xval[0] = m_tkrHists[ihist]->GetBinCenter(bin);
                    if(testAgainstFit) { 
                        yval = std::max(0.0, fitfunc->EvalPar(xval)); 
                    } else { 
                        yval = average[ihist];
                    }
                    bool deadBlock = (chipCount[ihist][(bin-1)/stripsPerBlock]==0);
                    bool enoughHits = yval >minHistCount || !checkCount;

                    // cutval based on the fitted value in the bin
                    if (testAgainstFit) {
                        cutval[0] = std::max(0., yval - nSig*sqrt(yval));
                        cutval[1] = std::max(hotFactor*yval, yval + nSig*sqrt(yval));
                    }

                    flagIt = false;
                    if (itest==0) { flagIt = deadBlock || xbin<=cutval[0]&& enoughHits;}
                    else          { flagIt = xbin>cutval[1]; }
                    
                    if (flagIt) {
                        testCount = true;
                        if (debugf) std::cout << "t,p,b,xb,cutval,f " << tower 
                            << " " << plane << " " << bin 
                            << " " << xbin << " " << cutval[itest] << " " 
                            << flagIt << std::endl;
                    }
                }
            }

            if (testCount) {
                if (inTower!=tower) {
                    if(!firstTower) listout << indent1.c_str() << "</tower>" << std::endl;
                    listout << indent1.c_str() << "<tower row=\"" << towerRow 
                        << "\" col=\"" << towerCol << "\"" << ">" << std::endl;
                }
                inTower = tower;
                firstTower = false;
                listout << indent2.c_str() << "<uniplane tray=\"" << tray << "\" which=\"";
                if (debugf) std::cout << "p,tr,tb" << plane << " " << tray << " " <<botTop << std::endl;
                if (botTop==0) listout << "bot" ;
                if (botTop==1) listout << "top" ;
                listout << "\" " << howBad;                               

                if (average[ihist]==0 && itest==0){
                    // do all-bad plane here
                    listout << " allBad=\"true\" />" << std::endl;                                
                } else {        
                    int first = -1; int last;
                    //bool inStripList = false;
                    bool isGap   = false;
                    float xbin;
                    int thisBin;
                    for (bin=1; bin<nbin+2; bin++) {
                        thisBin = bin;
                        if (bin==nbin+1) {
                            xbin = minHistCount+1;
                            xval[0] = m_tkrHists[ihist]->GetBinCenter(thisBin)+1;
                            yval = minHistCount+1;          
                        } else {
                            xbin = m_tkrHists[ihist]->GetBinContent(bin);
                            xval[0] = m_tkrHists[ihist]->GetBinCenter(bin);
                            if (testAgainstFit)  {
                                yval = std::max(0.0, fitfunc->EvalPar(xval));
                            } else {
                                yval = average[ihist];
                            }
                        }
                        bool enoughHits = yval>minHistCount || !checkCount;
                        int bin1 = (bin==nbin+1 ? nbin : bin);
                        bool deadBlock = (chipCount[ihist][(bin1-1)/stripsPerBlock]==0&&bin<nbin+1);

                        if (testAgainstFit) {
                            cutval[0] = std::max(0., yval - nSig*sqrt(yval));
                            cutval[1] = std::max(hotFactor*yval, yval + nSig*sqrt(yval));
                        }

                        flagIt = false;
                        if (itest==0) { flagIt = deadBlock || xbin<=cutval[0]&& enoughHits;}
                        else          { flagIt = xbin>cutval[1]; }
                        if (!flagIt && first>-1) isGap = true;
                   
                        if (flagIt) {
                            if(first==-1) {first = thisBin-1; last = first;}
                            else {if(!isGap) last++; }
                            if (newHist) {                                
                                listout << " >" << std::endl;
                                newHist = false;
                            }
                        }

                        if (isGap && first>-1) {
                            isGap = false;
                            if (last-first<10) {
                                if (badCount==0) {
                                    listout << indent3.c_str() <<  "<stripList strips= \" "; 
                                }
                                if (last-first==0) {
                                    listout << first << " ";
                                    if (debug) std::cout << first << " ";
                                    badCount++;
                                    totalBadCount++;
                                    if (badCount>0 && badCount%10==0) listout << std::endl  << "                    "; 
                                } else {
                                    int item;
                                    for (item=first; item<=last; ++item) {
                                        listout << item << " ";
                                        badCount++;
                                        totalBadCount++;
                                        if (badCount>0 && badCount%10==0) listout << std::endl  << "                    ";
                                    }
                                }
                            } else {
                                if (badCount) {
                                    listout << "\" />" << std::endl;
                                }
                                listout << indent3.c_str() << "<stripSpan first= \"" << first 
                                    << "\" last= \"" << last << "\" />" << std::endl;
                                badCount = 0;
                            }
                            first = -1;
                        }
                    }
                }   
                if (badCount ) listout << "\" />" << std::endl; 
                if (totalBadCount==0) listout << ">" ;
                if (testCount && !(average[ihist]==0&&itest==0)) listout << indent2.c_str() << "</uniplane>" << std::endl;
            }
            if (debugf) std::cout << "hist " << ihist << " done"<<std::endl;         
        }
        if (!firstTower) listout << indent1.c_str() << "</tower>" << std::endl;
        listout << "</badStrips>" << std::endl;
        listout.close();
    }

//Can't draw histos when running in compiled mode
#if defined(__CINT__)
    TCanvas* t;
    for (ihist=0; ihist<numHists; ++ihist){
        int pad = ihist%6 + 1;
        if(pad==1) {
            t = new TCanvas();
            gPad->Divide(2,3);
        }
        t->cd(pad);
        gPad->SetLogy(1);
        m_tkrHists[ihist]->Draw("Hist");
    } 

#endif
}
