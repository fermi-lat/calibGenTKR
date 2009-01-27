# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/calibGenTKR/SConscript,v 1.2 2008/11/04 14:54:33 ecephas Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: calibGenTKR-04-05-00
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




