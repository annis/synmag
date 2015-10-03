import numpy as np

"""Classes for describing DECam's varying system repsonse

"""

__author__ = ("Jim Annis <annis@fnal.gov> ")


class DECam(object):
    """Load system responses for all 62 ccds.
    For a given SED, predict magnitudes and colors.
    
    Public:
        mag(filter) 
        color(filter1, filter2) 
        magRel(filter, specialCCD = 30) 
    
    Private:
        loadSysResp () 
    """

    def __init__(self, 
            filterList = ["g","r","i","z","y"]
        ):
        """
        """
        self.filterList         = filterList
        masterDir = "/home/s1/annis/daedalean/synmag/sysresp/"
        self.masterDir          = masterDir
        self.sysResponseFileDir = masterDir
        self.sedDir             = masterDir + "../Stars/PICKLES/"
        self.outputDir          = masterDir 
        self.synMagFile         = self.outputDir + "synMags.txt"

        self.decam              = self.loadSysResp()


    def computeSynMagFile (self) :
        """
        """
        import os
        outputDir = self.outputDir
        synMagFile = self.synMagFile

        if (not os.path.exists(outputDir)) :
            os.mkdir(outputDir)

        picklesList = self.listOfPickles()
        picklesDict = self.picklesDict()
        filterList = self.filterList
        fd = open(synMagFile,"w")
        fd.write("# synthetic magnitudes of the pickles library\n")
        fd.write("# \t using the 2013 Wester DECam system response\n")
        fd.write("# \t and the Ting Li 1.3 airmass CTIO atmosphere\n")
        fd.write("# pickles-id, ")
        for i in range(1,63) :
            fd.write("ccd-{},".format(i))
        fd.write("filter\n")
        for star in picklesList :
            starnum = picklesDict[star]
            self.set_sed(star,"f_lambda", 10.)
            for filter in filterList :
                fd.write("{:3d} ".format(int(starnum)))
                mags = self.mag(filter)
                for ccd in range(1,63) :
                    fd.write("{:7.3f} ".format(mags[ccd-1]))
                filtnum = self.filterNumFromFilter(filter)
                fd.write("{:2d} \n".format(int(filtnum)))
        fd.close()

    def filterNumFromFilter(self, filter) :
        if filter == "u" : filtNum = 0
        elif filter == "g" : filtNum = 1
        elif filter == "r" : filtNum = 2
        elif filter == "i" : filtNum = 3
        elif filter == "z" : filtNum = 4
        elif filter == "y" : filtNum = 5
        return filtNum

    def set_sed(self,sedfile, sed_type, scale_lambda= 1.0) :
        """Input spectra can be:
            ergs/hz/sec/cm^2     = "f_nu"
            ergs/A/sec/cm^2      = "f_lambda"
            photons/hz/sec/cm^2  = "photon_nu"
            photons/A/sec/cm^2   = "photon_lambda"

        The Pickles library is in "f_lambda"
        One can scale the wavelengths to be in nm using scale_lambda.
            Pickles needs scale_lambda = 10.0
        """
        dir = self.sedDir
        sedfile = dir + sedfile
        self.sedFile = sedfile

        filterList         = self.filterList
        decam              = self.decam

        for filter in filterList :
            decam[filter].set_sed(sedfile, sed_type, scale_lambda)
    
    def colorRel(self, filter1, filter2, specialCCD = 36) :
        """ Calculate the color
        relative to a special CCD, by default #36.
        in all 62 CCDS
        """
        mags1 = self.mag(filter1)
        mags2 = self.mag(filter2)
        s_mag1 = mags1[specialCCD-1]
        s_mag2 = mags2[specialCCD-1]
        colors  = (mags1-mags2) - (s_mag1 - s_mag2)
        return colors

    def magRel(self, filter, specialCCD = 36) :
        """ Calculate the magnitude
        relative to a special CCD, by default #36.
        in all 62 CCDS
        """
        mags = self.mag(filter)
        s_mag = mags[specialCCD-1]
        mags = mags - s_mag
        return mags

    def color(self, filter1, filter2) :
        """ Calculate a color in all 62 ccds
        """
        mags1 = self.mag(filter1)
        mags2 = self.mag(filter2)
        color = mags1-mags2
        return color

    def mag(self, filter) :
        """ Calculate the magnitude in all 62 ccds
        """
        decam = self.decam

        mags = []
        for ccdnum in range(1,63) :
            decam[filter].set_ccd(ccdnum)
            mag = decam[filter].magnitude()
            mags.append(mag)
        mags = np.array(mags)
        return mags

    def loadSysResp (self) :
        """This loads 62 system response files into the object.
        The system response files are assumed to live in sysReponseFileDir
        and to have name in the format "DES-{}.txt".format(filter), with 
        one wavelength col and 62 ccd cols
        
        """
        import synMagnitude
        filterList = self.filterList
        decam = {}

        for filter in filterList :
            sys = self.sysResponseFileDir + "DES-{}.txt".format(filter)
            decam[filter] = synMagnitude.synMag(sys)

        return decam

    def listOfPickles(self) :
        """
        """
        pickles = [ "a0iii.sed", "a0i.sed", "a0iv.sed", "a0v.sed",
            "a2i.sed", "a2v.sed", "a3iii.sed", "a3v.sed", "a47iv.sed", 
            "a5iii.sed", "a5v.sed", "a7iii.sed", "a7v.sed", "b0i.sed", 
            "b0v.sed", "b12iii.sed", "b1i.sed", "b1v.sed", "b2ii.sed", 
            "b2iv.sed", "b3iii.sed", "b3i.sed", "b3v.sed", "b57v.sed", 
            "b5iii.sed", "b5ii.sed", "b5i.sed", "b6iv.sed", "b8i.sed", 
            "b8v.sed", "b9iii.sed", "b9v.sed", "f02iv.sed", "f0iii.sed",
            "f0ii.sed", "f0i.sed", "f0v.sed", "f2iii.sed", "f2ii.sed", 
            "f2v.sed", "f5iii.sed", "f5i.sed", "f5iv.sed", "f5v.sed", 
            "f6v.sed", "f8i.sed", "f8iv.sed", "f8v.sed", "g0iii.sed", 
            "g0i.sed", "g0iv.sed", "g0v.sed", "g2i.sed", "g2iv.sed", 
            "g2v.sed", "g5iii.sed", "g5ii.sed", "g5i.sed", "g5iv.sed", 
            "g5v.sed", "g8iii.sed", "g8i.sed", "g8iv.sed", "g8v.sed",
            "k01ii.sed", "k0iii.sed", "k0iv.sed", "k0v.sed", "k1iii.sed",
            "k1iv.sed", "k2iii.sed", "k2i.sed", "k2v.sed", "k34ii.sed", 
            "k3iii.sed", "k3i.sed", "k3iv.sed", "k3v.sed", "k4iii.sed", 
            "k4i.sed", "k4v.sed", "k5iii.sed", "k5v.sed", "k7v.sed", 
            "m0iii.sed", "m0v.sed", "m10iii.sed", "m1iii.sed", "m1v.sed", 
            "m2iii.sed", "m2i.sed", "m2p5v.sed", "m2v.sed", "m3iii.sed",
            "m3ii.sed", "m3v.sed", "m4iii.sed", "m4v.sed", "m5iii.sed",
            "m5v.sed", "m6iii.sed", "m6v.sed", "m7iii.sed", "m8iii.sed",
            "m9iii.sed", "o5v.sed", "o8iii.sed", "o9v.sed", "rf6v.sed",
            "rf8v.sed", "rg0v.sed", "rg5iii.sed", "rg5v.sed", "rk0iii.sed",
            "rk0v.sed", "rk1iii.sed", "rk2iii.sed", "rk3iii.sed", "rk4iii.sed",
            "rk5iii.sed", "wf5v.sed", "wf8v.sed", "wg0v.sed", "wg5iii.sed",
            "wg5v.sed", "wg8iii.sed", "wk0iii.sed", "wk1iii.sed", "wk2iii.sed",
            "wk3iii.sed", "wk4iii.sed", ]
        return pickles

    def picklesDict (self) :
        """ Look up table from pickles name to number,
        and vice versa
        """
        data =self.listOfPickles()
        pickles = dict()
        for i in range(0,len(data) ) :
            pickles[i] = data[i]
            pickles[data[i]] = i
        return pickles

# end namespace DECam
