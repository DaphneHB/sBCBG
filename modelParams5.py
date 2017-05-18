#!/apps/free/python/2.7.10/bin/python 

# defines the value of the parameters that will be used by testFullbG.py
# should be generated by sangoScript.py

params = {'LG14modelID':5,
          'nbCh':      2,
          'nbcpu':     6,
          'nbMSN': 2644.,
          'nbFSI':   53.,
          'nbSTN':    8.,
          'nbGPe':   25.,
          'nbGPi':   14.,
          'nbCSN': 3000.,
          'nbPTN':  100.,
          'nbCMPf':   9.,
          'GMSN':     4.5,
          'GFSI':     1.1,
          'GSTN':     1.35,
          'GGPe':     1.,
          'GGPi':     1.,
          'IeGPe':   13.,
          'IeGPi':   11.,
          'inDegCSNMSN': 100.,
          'inDegPTNMSN':   2.,
          'inDegCMPfMSN':  1.,
          'inDegFSIMSN':  30., # 30 : according to Humphries et al. 2010, 30-150 FSIs->MSN
          'inDegMSNMSN':  70., # 70 = 210/3 : according to Koos et al. 2004, cited by Humphries et al., 2010, on avg 3 synpase per MSN-MSN connection
          'inDegSTNMSN':   0.,
          'inDegGPeMSN':   0.,
          'inDegCSNFSI':  50.,
          'inDegPTNFSI':   17.,
          'inDegSTNFSI':   2.,
          'inDegGPeFSI':  40.,
          'inDegCMPfFSI':  9.,
          'inDegFSIFSI':  15., # 15 : according to Humphries et al., 2010, 13-63 FSIs->FSI
          'inDegPTNSTN':  25.,
          'inDegCMPfSTN':  9.,
          'inDegGPeSTN':  25.,
          'inDegCMPfGPe':  9.,
          'inDegSTNGPe':   8.,
          'inDegMSNGPe':2644.,
          'inDegGPeGPe':  25.,
          'inDegMSNGPi':2644.,
          'inDegSTNGPi':   8.,
          'inDegGPeGPi':  23.,
          'inDegCMPfGPi':  9.,
          }
