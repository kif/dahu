import math as m
import numpy as np
import sys

eps=3.66
el=float('1.6022e-19')
densSi=2.33
t=0.03
xp=[2.,3.,4.,5.,6.,8.,10.,15.,20.,30.,40.,50.,60.,80.,100.,150.,200.,300.,400.,500.] # mass energy-absorption coefficients
fp=[2669.,951.6,442.7,240.,143.9,63.13,32.89,9.794,4.076,1.164,0.4782,0.243,0.1434,0.06896,0.04513,0.03086,0.02905,0.02932,0.02968,0.02971]
# mass energy-absorption coefficients


def calcApe(E):
    return m.pow(10,(4.158-2.238*m.log10(E)-0.477*(m.log10(E)**2)+0.0789*(m.log10(E)**3)))
    
def calcApeEnergy(E):
    return np.interp(E,xp,fp)
    #return 24118.63*m.pow(E,-2.907)

def main(Int,E,atten='energy'):
    if atten=='photo':
        flux=eps*Int*float('1e-12')/(el*E*1000*(1-m.exp(-1*calcApe(E)*t*densSi)))
        print 'Flux=',flux,'photons/s, photoelectric absorption cross section'
    else:
        flux=eps*Int*float('1e-12')/(el*E*1000*(1-m.exp(-1*calcApeEnergy(E)*t*densSi)))
        print 'Flux=',flux,'photons/s, energy absorption cross section'
    return flux
   
		
if __name__ == "__main__":
   #print sys.argv
   if len(sys.argv)==4:
       main(float(sys.argv[1]),float(sys.argv[2]),str(sys.argv[3]))
   elif len(sys.argv)==3: 
       main(float(sys.argv[1]),float(sys.argv[2]),)
   else:
       print "---------"
       print "Calculates the beam intensity from photodiode reading"
       print "Usage: python calculate_energy.py photodiode_current(pA) Energy(keV)"  
       print "---------"
