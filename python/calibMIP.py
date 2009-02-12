
import sys, os, math, time, array, xml.dom.minidom

sys.path += [( os.path.join( os.environ.get("ROOTSYS"), "lib" ) )]
sys.path += [( os.path.join( os.environ.get("CALIBTKRUTILROOT"), "python" ) )]

import ROOT
import tkrUtils

ROOT.gSystem.Load("tkrPyRoot")
ROOT.gStyle.SetPalette(1)
#ROOT.gStyle.SetOptStat(0)

# get tag and version numbers
__tag__  = "$Name:  $"
__version__  = "$Revision: 1.4 $"
tagv = "%s:%s" % (__tag__.split()[1], __version__.split()[1])

params = { "minEntries":200.0, "minFracBadTot":0.08, "peakMIP":4.92, \
           "minGSigma":0.6, "maxGSigma":1.4, \
           "minLWidth":0.43, "maxLWidth":0.52, \
           "minFitGSigma":0.4, "maxFitGSigma":1.4, \
           "minFitLWidth":0.2, "maxFitLWidth":0.7, \
           "vRSigma":4.0, "vGFrac":0.78 }


def getDirName( fpath ):
  paths = fpath.split( '/' )
  dirName =  "/"
  for path in paths[:-1]:
    dirName = os.path.join( dirName, path )
  return dirName


class calibMIP:
  def __init__(self, jobOption, jname ):

    if not self.getJob( jobOption, jname ):
      sys.exit( "no job named %s found in %s" % (jname, jobOption) )

    self.timeStamps = []
    self.endTime = -1
    self.startTime = -1
    
    self.initHists()
    self.readTotParam()
    

  def getJob(self, jobOption, jname ):

    print "read job option: %s" % jobOption
    jobdir = getDirName( jobOption )
    print "job option dir: %s" % jobdir
    dom = xml.dom.minidom.parse( jobOption )
    topElm = dom.getElementsByTagName("jobOptions")[0]
    try:
      if jname == None: jname = topElm.getAttribute("default")
    except:
      sys.exit( "no job name is defined.")
    print "search for job: %s" % jname
      
    output = topElm.getElementsByTagName("output")[-1] # take last one
    self.outdir = output.getAttribute("dir")
    if not os.path.exists( self.outdir ):
      sys.exit( "invalid output dir: %s" % self.outdir )
    
    jobs = topElm.getElementsByTagName("jobOption")
    for job in jobs:
      name = str( job.getAttribute("name") )
      if name != jname: continue

      self.jName = jname
      self.jType = str( job.getAttribute("type") )
      self.jMode = str( job.getAttribute("mode") )
      self.timeStamp = time.strftime("%y%m%d-%H%M%S", time.gmtime() )
      fname = "ChargeScale-%s-%s.log" % (jname, self.timeStamp )
      lname = os.path.join( self.outdir, fname )
      print "log file:", lname
      self.logfile = open( lname, 'w' )
      self.logMessage( "calibMIP.py tag: %s" % tagv )
      self.logMessage( "tkrUtil.py tag: %s" % tkrUtils.tagv )
      self.getParams( topElm )
      self.getParams( job )

      elms = job.getElementsByTagName("totParam")
      self.totFiles = { "root":[], "xml":[] }
      for elm in elms:
        type = str( elm.getAttribute("type") )
        if not self.totFiles.has_key( type ):
          self.logMessage( "invalid totParam type: %s" % type )
          continue
        fnames = str(elm.getAttribute("files")).split() 
        for fname in fnames:
          if os.path.exists( fname ): self.totFiles[type].append( fname )
          else: self.logMessage( "invalid tot file path: %s"% fname )
      print "finished totParam elements"
         
      elms = job.getElementsByTagName("data")
      self.rnames = { "tree":[], "hist":[] }
      for elm in elms:
        type = str( elm.getAttribute("type") )
        if not self.rnames.has_key( type ):
          self.logMessage( "invalid data type: %s" %  type )
          continue
        fnames = str(elm.getAttribute("files")).split()
        for fname in fnames:
          fpath = os.path.join( jobdir, fname )
          if os.path.exists( fpath ): self.readFileList( type, fpath )
          elif os.path.exists( fname )  or fname[:5]=="root:":
            if fname[-5:] == ".root": self.rnames[type].append( fname )
            else: self.readFileList( type, fname )
          else: self.logMessage( "invalid data file path: %s" % fname )
      print "finished data elements"

    print "jobOption done"
    
    self.logMessage( "minEntries: %.0f" %  params["minEntries"] )
    self.logMessage( "minFracBadTot: %.2f" %  params["minFracBadTot"] )
    self.logMessage( "peakMIP: %.2f" %  params["peakMIP"] )
    names = [ "GSigma", "LWidth", "FitGSigma", "FitLWidth" ]
    for name in names:
      self.logMessage( "%s: %.2f - %.2f" \
                       % (name,params["min"+name],params["max"+name]) )
    self.logMessage( "RSigma, GFrac: %.2f, %.2f" \
                     % (params["vRSigma"], params["vGFrac"]) )
    
    return True


  def getParams(self, topElm):
    
    elms = topElm.getElementsByTagName("parameters")
    for elm in elms:
      for name in params.keys():
        if elm.hasAttribute(name):
          fval = float( elm.getAttribute(name) )
          if fval != params[name]:
            params[name] = fval
            self.logMessage( "new value for %s: %.2f" % (name,params[name]) )


    
  def logMessage(self, message ):
    print message
    self.logfile.write( message + "\n" )

    
  def initHists(self):
    self.hists = {}
    self.hists["LWidth"] = ROOT.TH1F( "LWidth", "Landau width", 100, 0.2, 0.7 )
    self.hists["GSigma"] = ROOT.TH1F( "GSigma", "Gauss sigma", 50, 0, 2.0 )
    self.hists["MPV"] = ROOT.TH1F( "MPV", "MIP MPV", 80, 3.0, 7.0 )
    
    self.hists["tot"] = ROOT.TH1F( "tot", "tot", 256, 0, 51.2 )
    self.hists["tower"] = ROOT.TH1F( "tower", "tower", 16, 0 , 16 )
    self.hists["unp"] = ROOT.TH1F( "unp", "unp", 36, 0, 36 )
    self.hists["strip"] = ROOT.TH1F( "strip", "strip", 1536, 0, 1536 )
    self.hists["charge"] = ROOT.TH1F( "charge", "charge", 200, 0, 20 )
    self.hists["rcharge"] = ROOT.TH1F( "rcharge", "raw charge", 200, 0, 20 )
    self.hists["ccharge"] = ROOT.TH1F( "ccharge", "corrected charge", 200, 0, 20 )
    self.hists["cfac"] = ROOT.TH1F( "cfac", "cfac", 40, 0.8, 1.0 )
    self.hists["p0TOT"] = ROOT.TH1F( "p0TOT", "TOT threshold", 150, 0.0, 3.0 )
    self.hists["p1TOT"] = ROOT.TH1F( "p1TOT", "TOT gain", 100, 0.0, 2.0 )
    self.hists["p2TOT"] = ROOT.TH1F( "p2TOT", "TOT quad", 100, 0.0, 0.2 )
    self.hists["chisqTOT"] = ROOT.TH1F( "chisqTOT", "TOT chisq", 100, 0, 5.0 )

    self.hists["chargeScale"] = ROOT.TH1F( "chargeScale", "charge scale", \
                                           100, 0, 2 )
    self.hists["entries"] = ROOT.TH1F( "entries", "Charge entries", \
                                       200, 0, 4000 )
    self.hists["fracBadTot"] = ROOT.TH1F( "fracBadTot", "fraction of bad TOT",\
                                       100, 0, 0.5 )
    self.hists["fitLWidth"] = ROOT.TH1F( "fitLWidth", "Landau width", \
                                         100, 0., 1.0 )
    self.hists["fitGSigma"] = ROOT.TH1F( "fitGSigma", "Gauss sigma", \
                                         50, 0, 2.0 )
    self.hists["fitProb"] = ROOT.TH1F( "fitProb", "Fut prob.", 100, 0, 1.0 )
    self.hists["fitChisqNdf"] = ROOT.TH1F( "fitChisqNdf", "Fut chi^2/NDF", \
                                           100, 0, 5.0 )

    ic = 1
    for key in ["rcharge","ccharge","charge"]:
      self.hists[key].SetLineWidth( 2 )
      self.hists[key].SetLineColor( ic )
      ic *= 2
    

  def readTotParam(self):
    self.totParams = [0]*tkrUtils.g_nTowers
    for tower in range(tkrUtils.g_nTowers):
      self.totParams[tower] = [0]*tkrUtils.g_nUniPlanes
      for unp in range(tkrUtils.g_nUniPlanes):
        self.totParams[tower][unp]= [0.0]*tkrUtils.g_nStrips

    for type in self.totFiles.keys():
      for fname in self.totFiles[type]:
        if type == "root":
          params = tkrUtils.readTotParamsFromRoot( fname )
        elif type == "xml":
          params = tkrUtils.readTotParamsFromXml( fname )
        for hw, pmap in params:
          tower = tkrUtils.g_Serials.index( hw )
          if tower<0 and tower>=tkrUtils.g_nTowers:
            self.logMessage( "invalid tower ID: %s in %s" % (hw, fname) )
            continue
          for unp in range(tkrUtils.g_nUniPlanes):
            for strip in range(tkrUtils.g_nStrips):
              self.totParams[tower][unp][strip] = ( pmap[0][unp][strip], \
                                                    pmap[1][unp][strip], \
                                                    pmap[2][unp][strip], \
                                                    pmap[3][unp][strip] )
              
    # look at distribution of parameters for bookkeeping purpose
    for tower in range(tkrUtils.g_nTowers):
      for unp in range(tkrUtils.g_nUniPlanes):
        for strip in range(tkrUtils.g_nStrips):
          (p0, p1, p2, chisq) = self.totParams[tower][unp][strip]
          #print tower, unp, strip, self.totParams[tower][unp][strip]
          self.hists["p0TOT"].Fill( p0 )
          self.hists["p1TOT"].Fill( p1 )
          self.hists["p2TOT"].Fill( p2 )
          self.hists["chisqTOT"].Fill( chisq )
          

  def readFileList(self, type, fname ):
    
    self.logMessage( "read file list: %s" % fname )
    file = open( fname )
    lines = file.readlines()
    file.close()

    for line in lines:
      sarray = line.split()
      if len(sarray) == 0: continue
      rname = sarray[0]
      if os.path.exists( rname ) or rname[:5]=="root:":
        self.rnames[type].append( rname )
      else: self.logMessage( "invalid data file path:", rname )
    self.logMessage( "total# of %s files: %d" \
                     % (type, len(self.rnames[type])) )
    

  def analyze(self, maxEntries=0 ):
    self.maxEntries = maxEntries
    self.entries = 0
    
    self.rcharge = [0]*tkrUtils.g_nTowers
    for tower in range(tkrUtils.g_nTowers):
      self.rcharge[tower] = [0]*tkrUtils.g_nUniPlanes
      for unp in range(tkrUtils.g_nUniPlanes):
        self.rcharge[tower][unp]= [0]*tkrUtils.g_nFE
        for fe in range(tkrUtils.g_nFE):
          self.rcharge[tower][unp][fe] = []
          self.rcharge[tower][unp][fe].append( [] )

    t0 = time.time()
    for rname in self.rnames["tree"]:
      self.analyzeTreeFile( rname )

      etime = time.time() - t0
      rate = self.entries /  etime
      self.logMessage( "%d processed in %.1f s, %.0f entries/s" \
                       % (self.entries, etime, rate ) )
      
      if self.maxEntries>0 and self.entries>=self.maxEntries: break
    if self.entries < 1E6: self.jMode = "hist"
    
    

  def analyzeTreeFile(self, rname ):
    self.logMessage( "analyze Tree root file: %s" % rname )
    rf = ROOT.TFile.Open( rname )
    if rf == None or not rf.IsOpen():
      self.logMessage( "root file access error, skip this file" )
      return

    # check fit parameter of chargeAll
    try:
      hist = rf.FindObjectAny( "chargeAll" )
      func = hist.GetFunction( "langau2" )
      LWidth = func.GetParameter(0)
      GSigma = func.GetParameter(3)
    except:
      self.logMessage( "error in fit function access for chargeAll" )
      return      
    self.hists["LWidth"].Fill( LWidth )
    self.hists["GSigma"].Fill( GSigma )
    if LWidth < params["minLWidth"] or LWidth > params["maxLWidth"]:
      self.logMessage( "LWidth %.2f is out of range: %.2f - %.2f" \
                       % (LWidth, params["minLWidth"], params["maxLWidth"]) )
      return
    if GSigma < params["minGSigma"] or GSigma > params["maxGSigma"]:
      self.logMessage( "GSigma %.2f is out of range: %.2f - %.2f" \
                       % (GSigma, params["minGSigma"], params["maxGSigma"]) )
      return
    MPV = func.GetParameter(1)
    self.hists["MPV"].Fill( MPV )

    tree = rf.FindObjectAny( "totInfo" )
    try: tree.GetEntries()
    except:
      self.logMessage( "no totInfo in this file" )
      return

    # get start time and end time
    treeTS = rf.FindObjectAny( "timeStamps" )
    endTime = -1
    firstRunId = -1
    lastRunId = -1
    for ien in range(treeTS.GetEntries()):
      treeTS.GetEntry( ien )
      if endTime < treeTS.endTime: endTime = treeTS.endTime
      if lastRunId < treeTS.lastRunId: lastRunId = treeTS.lastRunId
      if firstRunId < 0 or firstRunId > treeTS.firstRunId:
        firstRunId = treeTS.firstRunId
    date =  time.strftime("%Y/%m/%d, %H:%M", \
                          tkrUtils.getGMT( firstRunId ) )
    self.logMessage( "run ID: %d, %d, end time: %d, date: %s" \
                     % (firstRunId, lastRunId, endTime, date ) )
    self.timeStamps.append( (firstRunId, lastRunId, endTime) )
    if self.startTime < 0 or firstRunId < self.startTime:
      self.startTime = firstRunId
    if endTime > self.endTime: self.endTime = endTime
    treeTS.Reset()
    treeTS.Delete()
      
    for ien in range(tree.GetEntries()):
      tree.GetEntry( ien )
      rawTOT = tree.rawTOT
      tot = rawTOT % 256 * 0.2 + 1E-10 # convert to time
      self.hists["tot"].Fill( tot )
      id = rawTOT / 256
      strip = id % tkrUtils.g_nStrips
      id /= tkrUtils.g_nStrips
      unp = id % tkrUtils.g_nUniPlanes
      tower = id / tkrUtils.g_nUniPlanes
      fe = strip / tkrUtils.g_nChannels
      self.hists["tower"].Fill( tower )
      self.hists["unp"].Fill( unp )
      self. hists["strip"].Fill( strip )

      icharge = tree.charge
      charge = (icharge%100000) * 1.0E-3
      cfac = (icharge/100000) * 1.0E-3
      self.hists["charge"].Fill( charge )
      self.hists["cfac"].Fill( cfac )

      (p0, p1, p2, chisq) = self.totParams[tower][unp][strip]
      charge = cfac * ( p0 + (p1 + p2*tot) * tot  )
      if len( self.rcharge[tower][unp][fe][-1] ) == 100:
        self.rcharge[tower][unp][fe][-1] = array.array( 'f', self.rcharge[tower][unp][fe][-1] )
        self.rcharge[tower][unp][fe].append( [] )
      self.rcharge[tower][unp][fe][-1].append( charge )

      self.entries += 1
      if self.maxEntries>0 and self.entries>=self.maxEntries: break
    rf.Close()

  def fit(self):
    print self.jMode
    if self.jMode != "fit":
      self.saveHists()
      return

    self.logMessage( "start fit" )
    self.ffit = ROOT.defLangau( "langau", 0, 30 )
    self.ffit.SetParNames( "LWidth", "MP", "Area", "GSigma" )

    self.nerr = 0
    self.nwarn = 0
    self.hcharge = []
    self.hccharge = []
    self.chargeScale = []
    for tower in range(tkrUtils.g_nTowers):
      self.hcharge.append( [] )
      self.hccharge.append( [] )
      self.chargeScale.append( [] )
      for unp in range(tkrUtils.g_nUniPlanes):
        self.hcharge[tower].append( [] )
        self.hccharge[tower].append( [] )
        self.chargeScale[tower].append( [1.0]*tkrUtils.g_nFE )
        for fe in range(tkrUtils.g_nFE):
          lname = tkrUtils.g_layerNames[unp] 
          name = "chargeT%02dL%02dFE%02d" % (tower, unp, fe)
          title = "charge T%02d %s FE%02d" % (tower, lname, fe )
          hist = ROOT.TH1F( name, title, 200, 0, 20 )
          for ar in self.rcharge[tower][unp][fe]:
            for ic in range(len(ar)): hist.Fill( ar[ic] ) 
          self.hcharge[tower][unp].append( hist )
          (peak, peakErr, chisq, ndf) = self.fitTOT( tower, unp, fe )
          if peak <= 0.0:
            self.hccharge[tower][unp].append( None )
            self.chargeScale[tower][unp][fe] = (1.0,1.0,100.0,1)
            continue
          chargeScale = params["peakMIP"] / peak
          error = peakErr/peak * chargeScale
          self.chargeScale[tower][unp][fe] = (chargeScale,error,chisq,ndf)
          self.hists["chargeScale"].Fill( chargeScale )
          name = "cchargeT%02dL%02dFE%02d" % (tower, unp, fe)
          title = "corrected charge T%02d %s FE%02d" % (tower, lname, fe )
          hist = ROOT.TH1F( name, title, 200, 0, 20 )
          for ar in self.rcharge[tower][unp][fe]:
            for ic in range(len(ar)):
              #charge =  self.rcharge[tower][unp][fe][ic]
              charge =  ar[ic]
              hist.Fill( chargeScale * charge )
              self.hists["rcharge"].Fill( charge )
              self.hists["ccharge"].Fill( chargeScale * charge )
          self.hccharge[tower][unp].append( hist )
          self.rcharge[tower][unp][fe] = None

    mean = self.hists["chargeScale"].GetMean()
    rms = self.hists["chargeScale"].GetRMS()
    uplim = mean + 3*rms
    lowlim = mean - 3*rms
    for tower in range(tkrUtils.g_nTowers):
      for unp in range(tkrUtils.g_nUniPlanes):
        lname = tkrUtils.g_layerNames[unp] 
        for fe in range(tkrUtils.g_nFE):
          try:
            (chargeScale,error,chisq,ndf) = self.chargeScale[tower][unp][fe]
          except:
            print tower, unp, fe
            print len(self.chargeScale)
            print len(self.chargeScale[tower])
            print self.chargeScale[tower][unp]
          if chargeScale > uplim:
            self.nwarn +=1
            fename = "T%02d %s FE:%02d" % (tower, lname, fe)
            self.logMessage( "%s: chargeScale %.2f is greater than %.2f" \
                             % (fename,chargeScale,uplim) )
            self.chargeScale[tower][unp][fe] = (uplim,error,-chisq,-ndf)
          if chargeScale < lowlim:
            self.nwarn +=1
            fename = "T%02d %s FE:%02d" % (tower, lname, fe)
            self.logMessage( "%s: chargeScale %.2f is less than %.2f" \
                             % (fename,chargeScale,lowlim) )
            self.chargeScale[tower][unp][fe] = (lowlim,error,-chisq,-ndf)
          
    self.logMessage( "# of error: %d, # of warning: %d" \
                     % (self.nerr,self.nwarn) )

    func = ROOT.defLangau2( "langau2", 0, 30 )
    func.SetParNames( "LWidth", "MP", "Area", "GSigma", "RSigma", "GFrac" )
    keys = ["rcharge","charge","ccharge"]
    for key in keys:
      hist = self.hists[ key ]
      area = hist.Integral()
      mean = hist.GetMean()
      rms = hist.GetRMS()
      func.SetParLimits( 0, 0.0, rms )
      func.SetParLimits( 1, 0.0, mean*2 )
      func.SetParLimits( 2, 0.0, area*0.4 )
      func.SetParLimits( 3, 0.0, rms )
      func.SetRange( mean-1.25*rms, mean+2*rms )
      func.SetParameters(rms*0.5, mean*0.75, area*0.1, rms*0.4, \
                         params["vRSigma"], params["vGFrac"] )
      hist.Fit( "langau2", "RBQ" )
      self.logMessage( "*** %s fit results ***" % key )
      self.logMessage( "MPV: %.2f" % func.GetParameter(1) )
      self.logMessage( "LWidth: %.2f" % func.GetParameter(0) )
      self.logMessage( "GSigma: %.2f" % func.GetParameter(3) )
      self.logMessage( "RSigma: %.2f" % func.GetParameter(4) )
      self.logMessage( "GFrac: %.2f" % func.GetParameter(5) )
    
    self.saveHists()
    

  def fitTOT(self, tower, unp, fe):

    nwarn = 0
    chargeHist = self.hcharge[tower][unp][fe]
    entries = chargeHist.Integral()
    mean = chargeHist.GetMean()
    rms = chargeHist.GetRMS()
    self.hists["entries"].Fill( entries )
    lname = tkrUtils.g_layerNames[unp]
    fename = "T%02d %s FE:%02d" % (tower, lname, fe)

    if entries<params["minEntries"] or mean==0.0 or rms==0.0:
      alert = "%s, Entries %.0f, Mean: %.1f, RMS: %.1f skipped" \
              % (fename, entries, mean, rms)                
      self.logMessage( alert )
      self.nerr += 1
      return -1, -1, -1, -1

    binWidth = chargeHist.GetBinWidth( 2 )
    bin = int(mean*0.5/binWidth) + 1
    fracBadTot = chargeHist.Integral(1,bin) / entries
    self.hists["fracBadTot"].Fill( fracBadTot )

    lowLim = mean - 1.4 * rms
    if fracBadTot > params["minFracBadTot"] and lowLim < mean*0.5:
      lowLim = mean*0.5
      alert = "%s, large bad TOT fraction: %.3f > %.3f" \
              % (fename, fracBadTot, params["minFracBadTot"])
      self.logMessage( alert )
      nwarn += 1

    self.ffit.SetParLimits( 0, 0.1, 0.7 )
    self.ffit.SetParLimits( 1, 0.0, mean*2 )
    self.ffit.SetParLimits( 2, 0.0, entries*1.5 )
    self.ffit.SetParLimits( 3, rms*0.2, rms*1.5 )
    self.ffit.SetParameters( rms*0.2, mean*0.75, entries*0.1, rms*0.4 )
    self.ffit.SetLineColor(2)
    #self.ffit.SetRange( lowLim, mean+2.5*rms )
    #self.ffit.FixParameter( 4, m_RSigma )
    #self.ffit.FixParameter( 5, m_GFrac )
    chargeHist.Fit( "langau", "RBQ", "", lowLim, mean+2.5*rms )
    #canvas.Update()

    #0:width 1:peak 2:total entries 3:width(sigma)
    par = self.ffit.GetParameters()
    error = self.ffit.GetParErrors()
    self.hists["fitLWidth"].Fill( par[0] )
    self.hists["fitGSigma"].Fill( par[3] )
    self.hists["fitProb"].Fill( self.ffit.GetProb() )
    chisq = self.ffit.GetChisquare()
    ndf = self.ffit.GetNDF()
    self.hists["fitChisqNdf"].Fill( chisq/ndf )
    #self.towerHists[tower]["totPeak"].SetBinContent(unp+1,par[1])
    #self.towerHists[tower]["totPeak"].SetBinError(unp+1,error[1])
    if par[0] < params["minFitLWidth"] or par[0] > params["maxFitLWidth"]:
      alert = "%s, LWidth: %.2f is out of range %.2f - %.2f" \
              % (fename, par[0],params["minFitLWidth"],params["maxFitLWidth"])
      self.logMessage( alert )
      nwarn += 1
    if par[3] < params["minFitGSigma"] or par[3] > params["maxFitGSigma"]:
      alert = "%s, GSigma: %.2f is out of range %.2f - %.2f" \
              % (fename, par[3],params["minFitGSigma"],params["maxFitGSigma"])
      self.logMessage( alert )
      nwarn += 1

    if nwarn > 0: self.nwarn += 1
    return par[1], error[1], chisq, ndf


  def saveHists(self):

    fname = "ChargeScale-%s-%s.root" % (self.jName, self.timeStamp)
    rname = os.path.join( self.outdir, fname )
    self.logMessage( "save hist to %s" % rname )
    #rlist = ROOT.gDirectory.GetList()
    rf = ROOT.TFile( rname, "RECREATE" )
    for key in self.hists.keys():
      self.hists[key].Write()
    if self.jMode == "fit":
      for tower in range(tkrUtils.g_nTowers):
        tname = "T%02d" % tower
        ROOT.gDirectory.mkdir( tname )
        ROOT.gDirectory.cd( tname )
        for unp in range(tkrUtils.g_nUniPlanes):
          for fe in range(tkrUtils.g_nFE):
            self.hcharge[tower][unp][fe].Write()
            if self.hccharge[tower][unp][fe] != None:
              self.hccharge[tower][unp][fe].Write()
        ROOT.gDirectory.cd( ".." )
    #rlist.Write()

    # save time stampes into a tree
    endTime = array.array( 'f', [1] )
    firstRunId = array.array( 'L', [0] )
    lastRunId = array.array( 'L', [0] )
    tree = ROOT.TTree( "timeStamps", "timeStamps" )
    tree.Branch( "endTime", endTime, "endTime/f" )
    tree.Branch( "firstRunId", firstRunId, "firstRunId/i" )
    tree.Branch( "lastRunId", lastRunId, "lastRunId/i" )
    for firstRun, lastRun, eTime in self.timeStamps:
      #print firstRun, lastRun, eTime
      (firstRunId[0], lastRunId[0], endTime[0]) = (firstRun, lastRun, eTime)
      tree.Fill()
    tree.Write()
            
    rf.Close()

    startDate =  time.strftime("%Y/%m/%d, %H:%M", \
                               tkrUtils.getGMT( self.startTime ) )
    endDate = time.strftime("%Y/%m/%d, %H:%M", \
                            tkrUtils.getGMT( self.endTime ) )
    self.logMessage( "dates: %s - %s" %(startDate, endDate) )
    self.logMessage( "time duration: %d - %d" %(self.startTime,self.endTime) )
    deltaT = self.endTime - self.startTime
    if deltaT > 0:
      self.logMessage( "rate: %.1f entries/s" % (self.entries/deltaT) )

    if self.jMode == "fit":
      fname = "TkrChargeScale-%s-%s.root" % (self.jName, self.timeStamp)
      rname = os.path.join( self.outdir, fname )
      self.logMessage( "save chargeScale to %s" % rname )
      tkrUtils.saveChargeScaleRoot( rname, self.chargeScale )
    
    self.logMessage( "close log file: %s" % self.logfile.name )
    self.logfile.close()


if __name__ == '__main__':

  jname = None
  jobOption = os.path.join( os.environ.get('CALIBGENTKRROOT'), \
                            'python', "jobOptions.xml" )
  maxEntries = int( 1E5 )
  usage = "usage: python calibMIP [job name] -j [job option file] -m [max entries] -b"

  for iarg in range(1,len(sys.argv)):
    if sys.argv[iarg] == "-n": maxEntries = -1
    elif sys.argv[iarg] == "-j": jonOption = "TBD"
    elif sys.argv[iarg] == "-b": ROOT.gROOT.SetBatch()
    elif sys.argv[iarg] == "-h": sys.exit( usage )
    elif maxEntries < 0:
      maxEntries = int( float( sys.argv[iarg] ) )
      if maxEntries < 0: maxEntries = 0
    elif jobOption == "TBD": jobOption = sys.argv[iarg]
    elif jname == None: jname = sys.argv[iarg]
    else: sys.exit( usage )

  calib = calibMIP( jobOption, jname )
  calib.analyze( maxEntries )
  calib.fit()

  #calib.hchargeAll.SetLineColor( 4 )
  #calib.hchargeAll.Draw("same")
