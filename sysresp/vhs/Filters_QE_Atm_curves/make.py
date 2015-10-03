import numpy as np
from scipy.interpolate import interp1d

def them() :
    print "make vhs filters"
    qe = getQE()
    water = watervapor()
    tele = telescope()
    jfilter = jband()
    hfilter = hband()
    kfilter = kband()

    qe_interp = interp1d(qe[0], qe[1], bounds_error=False, fill_value=0.0)
    water_interp = interp1d(water[0], water[1], bounds_error=False, fill_value=0.0)
    tele_interp = interp1d(tele[0], tele[1], bounds_error=False, fill_value=0.0)

    j = jfilter[1]*qe_interp(jfilter[0])*water_interp(jfilter[0])*tele_interp(jfilter[0])
    h = hfilter[1]*qe_interp(hfilter[0])*water_interp(hfilter[0])*tele_interp(hfilter[0])
    k = kfilter[1]*qe_interp(kfilter[0])*water_interp(kfilter[0])*tele_interp(kfilter[0])

    np.savetxt("VHS_J.txt", np.array([jfilter[0], j]).T, "%.1f %.5f")
    np.savetxt("VHS_H.txt", np.array([hfilter[0], h]).T, "%.1f %.5f")
    np.savetxt("VHS_K.txt", np.array([kfilter[0], k]).T, "%.1f %.5f")


def getQE() :
    file = "qe.tab"
    wave,trans = np.genfromtxt(file, unpack=True)
    trans = trans/100.
    return wave, trans

def watervapor () :
    # 3 mm water, airmass = 2.0
    #file = "trans_30_20.dat"
    #wave,trans = np.genfromtxt(file, unpack=True,skiprows=2)

    # 1 mm water, airmass = 1.0
    file = "trans_10_10.dat"
    wave,trans = np.genfromtxt(file, unpack=True,skiprows=2)

    wave = wave*1000. ;# microns to nm
    return wave, trans

def telescope() :
    file = "VISTA_M1_Reflectivity_forETC_2009-Sep12.txt"
    wave,m1 = np.genfromtxt(file, unpack=True)
    file = "VISTA_M2_Reflectivity_forETC_2007-Jun19.txt"
    wave,m2 = np.genfromtxt(file, unpack=True)

    telescope = m1*m2
    telescope = telescope/100./100.
    return wave, telescope

def jband() :
    file="VISTA_Filters_at80K_forETC_J.dat"
    wave,filter = np.genfromtxt(file, unpack=True)
    filter = filter/100.
    return wave, filter
def hband() :
    file="VISTA_Filters_at80K_forETC_H.dat"
    wave,filter = np.genfromtxt(file, unpack=True)
    filter = filter/100.
    return wave, filter
def kband() :
    file="VISTA_Filters_at80K_forETC_Ks.dat"
    wave,filter = np.genfromtxt(file, unpack=True)
    filter = filter/100.
    return wave, filter

