# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/calibGenTKR/SConscript,v 1.7 2009/02/10 17:30:36 glastrm Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: calibGenTKR-04-08-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()


progEnv.Tool('addLibrary', library = baseEnv['rootGuiLibs'])
progEnv.Tool('addLibrary', library = baseEnv['rootLibs'])
progEnv.Tool('facilitiesLib')
progEnv.Tool('commonRootDataLib')
progEnv.Tool('xmlBaseLib')
progEnv.Tool('calibUtilLib')
progEnv.Tool('calibRootDataLib')
progEnv.Tool('TkrUtilLib')
progEnv.Tool('digiRootDataLib')
progEnv.Tool('reconRootDataLib')
progEnv.Tool('calibTkrUtilLib')

doBadStripsCalib = progEnv.Program('doBadStripsCalib', listFiles(['src/badStripsCalib/*.cxx']))
doMuonCalibTot = progEnv.Program('doMuonCalibTot',listFiles(['src/muonCalibTot/*.cxx']))
doXml2root = progEnv.Program('doXml2root',listFiles(['src/xml2root/*.cxx']))

progEnv.Tool('registerObjects', package = 'calibGenTKR', binaries = [doMuonCalibTot,doXml2root,doBadStripsCalib])





