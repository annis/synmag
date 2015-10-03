import numpy as np
import synMagnitude
import DECam

def computeSynMagFile (outfile="../sysresp/symMags-sdss.txt") :
        """
        """
        import os
        synMagFile = outfile

        sm=dict()
        sm[0]=synMagnitude.synMag("../sysresp/SDSS-u.txt")
        sm[1]=synMagnitude.synMag("../sysresp/SDSS-g.txt")
        sm[2]=synMagnitude.synMag("../sysresp/SDSS-r.txt")
        sm[3]=synMagnitude.synMag("../sysresp/SDSS-i.txt")
        sm[4]=synMagnitude.synMag("../sysresp/SDSS-z.txt")
        de = DECam.DECam()
        
        picklesList = de.listOfPickles()
        picklesDict = de.picklesDict()
        filterList = [0,1,2,3,4]
        print "opening   ", synMagFile
        fd = open(synMagFile,"w")
        fd.write("# synthetic magnitudes of the pickles library\n")
        fd.write("# pickles-id, ")
        fd.write("ccd-{},".format(1))
        fd.write("filter\n")
        for star in picklesList :
            print star, 
            starnum = picklesDict[star]
            for filter in filterList :
                sm[filter].set_sed("../Stars/PICKLES/" + star,"f_lambda", 10.)
                fd.write("{:3d} ".format(int(starnum)))
                mag = sm[filter].magnitude()
                print mag, " ",
                fd.write("{:7.3f} ".format(mag))
                fd.write("{:2d} \n".format(filter))
            print ""
        print "closing   ", synMagFile
        fd.close()
