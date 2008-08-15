# -*- python -*-
# $Header$
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: calibGenTKR-04-05-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('calibGenTKRLib', depsOnly = 1)

progEnv.Tool('calibGenTKRLib')

doBadStripsCalib = progEnv.Program('doBadStripsCalib', listFiles(['src/badStripsCalib/*.cxx']))
doMuonCalibTot = progEnv.Program('doMuonCalibTot',listFiles(['src/muonCalibTot/*.cxx']))
doXml2root = progEnv.Program('doXml2root',listFiles(['src/xml2root/*.cxx']))

progEnv.Tool('registerObjects', package = 'calibGenTKR', binaries = [doBadStripsCalib,doMuonCalibTot,doXml2root])

