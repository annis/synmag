import numpy as np
import scipy as sp
import scipy.interpolate

"""Classes for describing a system response.

   *** Initializing a sysResp runs code to generate systerm response curves

"""

__author__ = ("Jim Annis <annis@fnal.gov> ")

class sysResp(object):
    """system response
    aerosolData can be 1,2,3, and chooses which of the aerosol data fits
    of Gutierrez-Moreno et al 196 to use.
    A choice of 0 gives no change over the palomar atmopshere.

    aerosols are currently disabled.

    Run sysResp.go() to generate five filter system response files
    DES-g.txt, DES-r.txt, DES-i.txt, DES-z.txt, DES-y.txt
    """

    def __init__(self, 
        ccdpos = "ccdPos-sky-v6.par",
        insrespDir = "filter_curves/",
        airmass = 1.3,
        aerosolData = 0
        ):
        """Create a sysResponse"""

        masterDir = "/home/s1/annis/daedalean/synmag/sysresp/"
        self.masterDir = masterDir 
        self.ccdPosFile = masterDir + ccdpos
        self.insresponse = masterDir + insrespDir
        self.outDir = masterDir

        self.airmass = airmass
        self.altitude = 2200.0
        # ctio aerosol coefficients from Gutierrez-Moreno et al 1986
        A_h = [0.0, 0.05, 0.018, 0.013, 0.017]
        alpha = [0.0, 1.0, 1.1, 1.2, 1.8]
        self.A_h   = A_h[aerosolData]
        self.alpha = alpha[aerosolData]

        self.filterList = [1,2,3,4,5]
        self.ccdData =self.getCCDPositions() 


#=================================================================================
#
# Main Routines
#
#=================================================================================
    def go (self) :
        """
        for the five des filters, construct and save the system response
        """
        filterList = self.filterList
        atmoData = self.make_atmosphere()
        outDir = self.outDir

        trans = dict()
        for filter in filterList :
            fil = self.nameFromFilterNumber(filter)
            filename = outDir + "DES-{}.txt".format(fil)
            self.sysRespHeader(filename, fil) 

            for ccd in range(1,63) :
                wave, transSpline = self.get_instrument_response(filter, ccd) 
                wave, trans[ccd] = self.make_sys_ccd(atmoData[0],atmoData[1],wave,transSpline)
            print "\t writing {}".format(filename)
            self.sysRespWrite (filename, wave, trans) 


    def make_sys_ccd(self, atWaves, atTransSpline, ccdWaves, ccdTransSpline) :
        """
        Combine the atmospheric and instrumental response
        on a 1nm grid from 300 to 1100 nm range

        This assumes that both atmo and inst response are splines
        and their respective waves only cover the range of the spline
        """

        waves = np.arange(300,1101)
        trans = np.zeros(waves.size)

        ix = np.nonzero( (waves >= atWaves[0]) & (waves <= atWaves[-1]) )
        atmo = atTransSpline(waves[ix])
        trans[ix] = atmo

        ix = np.nonzero( (waves >= ccdWaves[0]) & (waves <= ccdWaves[-1]) )
        ccd = ccdTransSpline(waves[ix])
        trans[ix] = trans[ix]*ccd

        return waves, trans


#=================================================================================
#
# Support  Routines
#
#=================================================================================
    def get_instrument_response(self, filter, ccd) :
        """
        A routine to deal with the way William gave us instrument
        response in 2013
        """
        maxRadius = 1.1 ;# degees
        ccdData = self.ccdData 
        ccdnum, x,y = ccdData[0], ccdData[1], ccdData[2]

        insDir = self.insresponse 
        ins = dict()
        ins["g"] = insDir + "g_003.dat"
        ins["r"] = insDir + "r_005.dat"
        ins["i"] = insDir + "i_003.dat"
        ins["z"] = insDir + "z_003.dat"
        ins["y"] = insDir + "y_003.dat"

        fil = self.nameFromFilterNumber(filter)
        data = self.get_ins_filter(ins[fil])
        # t = average over all non-exluded amplifiers
        # tr1 = averaged only for the inner two CCDs (<10% RMax)
        # tr2 = for 10<Rmax<30%
        # tr3 = for 30<Rmax<60%
        # tr4 = for Rmax>60%
        ix = np.nonzero(ccdnum == ccd)
        radius = np.sqrt( x[ix]**2 + y[ix]**2 )/maxRadius
        wave = data[0]
        if radius <= 0.1 : trans = data[2]
        elif (radius >  0.1) & (radius <= 0.3) : trans = data[3]
        elif (radius >  0.3) & (radius <= 0.6) : trans = data[4]
        elif (radius >  0.6) : trans = data[5]
        return wave, trans

    def get_ins_filter(self, file) :
        wave, t, tr1, tr2, tr3, tr4 = np.genfromtxt(file, unpack=True)
        t = t/1000.
        tr1 = tr1/1000.
        tr2 = tr2/1000.
        tr3 = tr3/1000.
        tr4 = tr4/1000.
        st = sp.interpolate.InterpolatedUnivariateSpline(wave, t)
        str1 = sp.interpolate.InterpolatedUnivariateSpline(wave, tr1)
        str2 = sp.interpolate.InterpolatedUnivariateSpline(wave, tr2)
        str3 = sp.interpolate.InterpolatedUnivariateSpline(wave, tr3)
        str4 = sp.interpolate.InterpolatedUnivariateSpline(wave, tr4)
        return wave, st, str1, str2, str3, str4
        
    def make_atmosphere(self) :
        """
        Ting Li's atmosphere
        """
        file = self.masterDir + "atmo/ctio.txt"
        airmass = self.airmass
        altitude = self.altitude
        A_h   = self.A_h   
        alpha = self.alpha 

        # ctio.txt gives wave (nm) and T
        waves = []
        trans = []
        waves, trans = np.genfromtxt(file,unpack=True,comments="#")

        # as per Gunn sdss mailing list: 
        #   deals with Rayleigh. And Ozone, sort of.
        # transMag = transMag*airmass*np.e**((1700-altitude)/7000)
        # as per Gutierrez-Moreno et al 1982
        #       with coefficients from Gutierrez-Moreno et al 1986
        # wave_pivot = waves[29] ;# pivot around 400 nm
        # aerosols = (1.086 *A_h * (waves/wave_pivot)**-alpha) * airmass
        # transMag = transMag + aerosols
        # trans = 10**(-0.4 * transMag)

        atmosphere = sp.interpolate.InterpolatedUnivariateSpline(waves, trans)
        return waves, atmosphere

    def nameFromFilterNumber (self, filter) :
        if filter == 0 : fil = "u"
        elif filter == 1 : fil = "g"
        elif filter == 2 : fil = "r"
        elif filter == 3 : fil = "i"
        elif filter == 4 : fil = "z"
        elif filter == 5 : fil = "y"
        return fil

#=================================================================================
#
# Read data
#
#=================================================================================
    def getCCDPositions(self) :
        """ return positions in cm
            ccdno,x,y
        """
        ccdpos = self.ccdPosFile 

        ccdno = []
        x = []
        y = []
        file = open(ccdpos,"r")
        for line in file:
            if (not line.split())  : continue
            if (line.split()[0] != "CCD") : continue
            ccdno.append(int( line.split()[8] ))
            x.append(float( line.split()[4] ))
            y.append(float( line.split()[5] ))

        ccdno = np.array(ccdno)
        # x,y in arcminutes
        x = np.array(x)/60.
        y = np.array(y)/60.
        # if x,y wanted in centimeters: 
        # x = x*18.181818
        # y = y*18.181818
        return [ccdno, x,y]

#=================================================================================
#
# Write data
#
#=================================================================================
    def sysRespHeader (self, filename, filterName) :
        import datetime
        airmass  = self.airmass
        altitude = self.altitude
        A_h      = self.A_h   
        alpha    = self.alpha 
        insfile = "William Wester's 2013 DECam System (instrument) response curves"
        insfile2 = "http://home.fnal.gov/~wester/work/"
        atfile = "Ting Li's airmass=1.3 CTIO atmosphere"
        atfile2 = "uvspec_afglus_pressure780_airmass1.3_asea1_avul1_pw03_tau0.03.out.txt"

        now = datetime.datetime.now()

        f = open(filename,"w")
        f.write("# \n")
        f.write("# DES System response, filter {}\n".format(filterName))
        f.write("# \n")
        f.write("# \tJim Annis here and now      {}\n".format(
            now.strftime("%Y-%m-%d %H:%M")))
        f.write("# \n")
        f.write("# Wester's file has responses in 4 annuli.\n")
        f.write("# The instrument response in the annuli in which the ccd center lies\n")
        f.write("# was taken as the instrument response for that ccd.\n")
        f.write("# \n")
        f.write("# instrument transmission:  \n")
        f.write("#       {}\n".format(insfile))
        f.write("#       {}\n".format(insfile2))
        f.write("# atmosphere transmission:  \n")
        f.write("#       {}\n".format(atfile))
        f.write("#       {}\n".format(atfile2))
        f.write("# \n")
        f.write("# airmass               {} \n".format(airmass))
        f.write("# aerosols, A_h, alpha  {} {} (but Ting's atm includes aerosols)\n".format(A_h, alpha))
        f.write("# altitude              {}\n".format(altitude))
        f.write("# \n")
        f.write("# wavelength(nm) Trans(ccd 1) Trans(ccd 2)... Trans(ccd 62)\n")
        f.write("# \n")
        f.close()

    def sysRespWrite (self, filename, waves, trans) :
        fd = open(filename,"a")
        for i in range(0,waves.size) :
            fd.write("{:4.0f} ".format(waves[i]))
            for ccd in range(1,63) :
                fd.write("{:7.4f} ".format(trans[ccd][i]))
            fd.write("\n")
        fd.close()

# end sysresponse namespace
