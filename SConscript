# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/calibGenTKR/SConscript,v 1.12 2012/05/09 19:07:45 jrb Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: calibGenTKR-04-08-03
Import('baseEnv')
Import('listFiles')
Import('packages')
if baseEnv['PLATFORM'] != "win32":
  progEnv = baseEnv.Clone()

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

  doBadStripsCalib = progEnv.Program('doBadStripsCalib',
                                     listFiles(['src/badStripsCalib/*.cxx']))
  doMuonCalibTot = progEnv.Program('doMuonCalibTot',
                                   listFiles(['src/muonCalibTot/*.cxx']))
  doXml2root = progEnv.Program('doXml2root',listFiles(['src/xml2root/*.cxx']))

  progEnv.Tool('registerTargets', package = 'calibGenTKR',
               binaryCxts = [[doMuonCalibTot,progEnv], [doXml2root,progEnv],
                             [doBadStripsCalib,progEnv]])





