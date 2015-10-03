import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy  import integrate

"""Classes for describing a synthetic magnitude.

    synMag loads a single system response, and a single SED

"""

__author__ = ("Jim Annis <annis@fnal.gov> ")

class synMag(object):
    """system response convolved with object flux
    The integration is via simpsons rule: scipy.integrate.simps

    wave and flux are numpy arrays
    Wavelength is stored as nm, and conversions for given files
    are made possible
    
    The SED flux can be input as
            ergs/hz/sec/cm^2     = "f_nu"
            ergs/A/sec/cm^2      = "f_lambda"
            photons/hz/sec/cm^2  = "photon_nu"
            photons/A/sec/cm^2   = "photon_lambda"
    Having  said that, internally we'll work in photons/sec/cm^2/Hz

    Here we also adopt the notion that what is fundamental isn't
    ergs/sec/cm^2/Hz or ergs/sec/cm^2/A but ergs/sec/cm^2.
    This means that we plot 
        nu F_nu vs log(nu) or lambda F_lambda vs log(lambda)
    For photons,
        nu photons vs log(nu) or lambda photons vs log(lambda)
    Think of this as per fractional (lambda or nu ) interval
    So there are four additional plots:
        ergs/sec/cm^2        = "energy" =  lambda f_lambda vs log(lambda) 
        photons/sec/cm^2     = "photon" = lambda photon_lambda vs log(lambda)

    Public:
        set_ccd         ( ccdnum)
        set_sed         ( sedfile, sed_type, scale_lambda=1.0, scale_flux=1.0) 
        set_blackbody_sed ( temperature_in_kelvin, scale_flux=1.0 )
        set_sed_v       (sed_wave, sed_flux, sed_type) 
        magnitude       () 
        logFlux         () 
        interpolateSed  (wavelengths, sed_waves, sed_flux) 
        redshiftSed     (zed, waves, flux, sedType) 
        int_flambda      ()      ;# integrate over all f_lambda
        plot            (figure, plotSedType) 

   Private:
        loadSystemResponse() 

    Example:
        sm = synMagnitude.synMag("/home/s1/annis/daedalean/synmag/sysresp/DES-g.txt")
        sm.set_ccd(36)
        sm.set_sed("Stars/PICKLES/k3iii.sed", "f_lambda", 10.)
            # Pickles library is in f_lambda and needs 10. wavelength scaling
        sm.plot(figure,"photon_lambda")
        g_mag = sm.magnitude()

    There is a blackbody function:
        sm.set_blackbody_sed(5777.)
        sm.plot(plt,"f_nu")

    Sys repsonse files
        /home/s1/annis/daedalean/synmag/sysrep
    """

    def __init__(self, 
            sysResponseFile = ""
        ):
        """"""
        import numpy as np
        # cgs
        self.h = 6.6261*10**-27
        self.c = 2.99792458*10**10
        self.hc = 1.9864847852e-16

        self.sysResponseFile = sysResponseFile
        if "SDSS" in sysResponseFile :
            self.camera          = "sdss"
        if "DES" in sysResponseFile :
            self.camera          = "des"
        waves, trans         = self.loadSystemResponse()
        self.waves           = waves
        self.trans           = trans
        if "SDSS" in sysResponseFile :
            self.set_ccd(1)

    def set_ccd(self, ccd) :
        """ set a ccd number and hence a system response
        Side effects are to compute phi and lambda_e
        """
        self.qe              = self.trans[ccd]
        self.phi             = self.Phi()
        self.lambda_e        = self.Lambda_e()

    def Phi (self) :
        """ compute the filter times the wavelength
        """
        qe = self.qe
        waves = self.waves
        norm = scipy.integrate.simps(qe*waves, waves)
        phi             = qe*waves/norm
        return phi

    def Lambda_e (self) :
        """ compute the system responses effective wavelength
        """
        waves = self.waves
        phi  = self.phi
        ln_waves = np.log(waves)
        ln_lambda_e = scipy.integrate.simps(ln_waves*phi, waves)
        lambda_e = np.e**(ln_lambda_e)
        return lambda_e
#
# the spectrum computed is Int B_nu dnu 
#     recall that Int B_nu = du = sigma T^4/pi
#
    def set_blackbody_sed (self, temperature, scale_flux = 1.0) :
        """Input blackbody will be in f_nu,
            ergs/hz/sec/cm^2     = "f_nu"
        Internally we'll work in photons/Hz/sec/cm^2 
        """
        # sigma_b =  5.67e-5 erg cm-2 K-4 s-1
        c = 3.00e10 ;# cm/s
        h = 6.63e-27 ;# erg-s
        k = 1.38e-16 ;# erg/K
        const1 = 2*h/c**2
        const2 = h/k
        sed_type = "f_nu"

        sed_flux = []
        sed_waves = np.arange(200,5001,0.1)
        for wave in sed_waves :
            nu = c/(wave*1e-7) ;# nm -> cm
            flux = const1 * nu**3 / (np.e**(const2 * nu/temperature) - 1 )
        #    if wave == 500: 
        #        print const1*nu**3, const2*nu/temperature, flux

            # this one keeps interesting NS scaling
            # L in units of 10^{40} ergs/sec, d in units of 100 Mpc, T in units of 1000 K
            # so the above scaling is
            #T=2.255; d= 1.0; L = 1.38
            #flux = flux * 4.4e-22  * L/d**2/T**4
            flux = flux * scale_flux

            sed_flux.append(flux)
            
        sed_flux = np.array(sed_flux)
        self.set_sed_v ( sed_waves, sed_flux, sed_type) 

    def set_sed (self, sedfile, sed_type, scale_lambda=1.0, scale_flux=1.0) :
        """Input spectra can be:
            ergs/hz/sec/cm^2     = "f_nu"
            ergs/A/sec/cm^2      = "f_lambda"
            photons/hz/sec/cm^2  = "photon_nu"
            photons/A/sec/cm^2   = "photon_lambda"
        Internally we'll work in photons/Hz/sec/cm^2, check self.sed_type
        """
        sed_waves, sed_flux = np.loadtxt(
            sedfile, comments='#',unpack=True, usecols=(0,1))
        sed_waves = sed_waves/scale_lambda
        sed_flux = sed_flux*scale_flux 

        self.set_sed_v ( sed_waves, sed_flux, sed_type) 

    def set_sed_v (self, sed_waves, sed_flux, sed_type, test=0) :
        """Input spectra can be:
            ergs/hz/sec/cm^2     = "f_nu"
            ergs/A/sec/cm^2      = "f_lambda"
            photons/hz/sec/cm^2  = "photon_nu"
            photons/A/sec/cm^2   = "photon_lambda"
        Internally we'll work in photons/Hz/sec/cm^2, check self.sed_type

        We'll also interpolate to the wavelengths of the system response
        """
        h = self.h
        c = self.c
        hc = self.hc
        # what we'll do is to convert everything into f_nu,
        # then at the end, we'll choose to work in photons/Hz
        waves_cm = sed_waves.copy()*10**-7  #; nm->cm
        if (sed_type == "f_nu") :
            sed_fnu = sed_flux                  ;# f_nu -> f_nu
        elif (sed_type == "f_lambda") :
            sed_fnu = sed_flux*waves_cm**2/c    ;# f_lambda -> f_nu
        elif (sed_type == "photon_nu") :
            sed_fnu = sed_flux*hc/waves_cm**1  ;# photon_nu -> f_nu
        elif (sed_type == "photon_lambda") :
            sed_fnu = sed_flux*h*waves_cm**1  ;# photon_lambda -> f_nu
        else :
            raise NameError('sed type not known: {}'.format(sed_type))

        # now we'll convert to photons/s/cm^2/Hz
        sed_p_nu = sed_fnu*waves_cm/hc
        
        # and we'll interpolate to the positions of the phi
        sed_p_nu_interp = self.interpolateSed(self.waves, sed_waves, sed_p_nu)

        # JTA testing
        #test = 1
        if test:
            print "Testing... "
            print "\t Blanco gets ~0.01 photons/s/cm^2 at 20th mag",
            print "which is 3.631e-28 ergs/sec/cm^2/Hz"
            # convert the spectra to have a 20th magnitude
            EperP = hc/(self.waves*10**-7.)
            fnu = EperP*sed_p_nu_interp 
            fnu_mean = fnu.mean()
            #fnu_mean = 1.0
            print "fnu_mean init = ", fnu_mean
            fnu = (fnu/fnu_mean) * 3.631*10**-28 ;# -> 20th mag
            #fnu = (fnu/fnu_mean) 
            print "fnu_mean final = ", fnu.mean()
            sed_p_nu_interp = fnu/EperP
            print "p_nu_mean  = {:.3e}   mean E per photon= {:.3e}".format(
                 sed_p_nu.mean(), EperP.mean())


        self.sed = sed_p_nu_interp
        self.sed_raw_waves = sed_waves
        self.sed_raw = sed_p_nu
        self.sed_type = "photon_nu"

    # integral of flambda from wave min to wave max
    # useful in normalizing spectra
    def int_flambda (self) :
        h = self.h
        c = self.c
        p_nu = self.sed_raw
        waves = self.sed_raw_waves

        waves_cm = waves.copy()*10**-7 #; nm->cm
        flambda = p_nu * h *c**2 / waves_cm**3
        integral = scipy.integrate.simps(flambda, waves_cm) 
        return integral


    def magnitude ( self) :
        """ Death to -2.5
        """
        logFlux = self.logFlux( )
        magnitude = -2.5*logFlux 
        return magnitude

    def logFlux ( self ) :
        """ Convert the measured number of photons to an energy flux
        ergs/sec/cm^2/Hz

        lambda_e is in what ever the units of the sys response is in,
        but the conversion here assumes it is in nm.
        """
        hc = self.hc   ;# speed of light  * plank constant
        lambda_e = self.lambda_e 
        logN_photons = self.logN_photons_mean()
        # these are the photons + conversion to energy+ units to cm from nm,
        # and the 48.60 AB zeropoint with a -2.5log(hc) in it.
        logEnergyFlux = logN_photons + np.log10(hc) \
            - np.log10(lambda_e) - np.log10(100./10**9.) \
            + 48.60/2.5
        return logEnergyFlux

    def logN_photons_mean ( self ) :
        """Convolve an SED, "sed", with a system reponse curve
        and return the number of detected photons, in the monochromatic
        AB sense.

        The way to think about this is that the convolution is
        aiming to get the weighted average photon flux.

        sed is in photons/s/cm^2/Hz
        thus this is an AB magnitude
            if the spectra was in /A, it would be a ST magnitude
        """

        phi = self.phi 
        sed = self.sed
        waves = self.waves
        photons = scipy.integrate.simps(sed*phi, waves)
        logN_photons = np.log10(photons)
        return logN_photons

    def logN_photons ( self ) :
        """Convolve an SED, "sed", with a system reponse curve
        and return the number of detected photons.

        sed is in photons/s/cm^2/Hz
        """

        # load the data
        c = self.c
        sed = self.sed
        phi = self.phi 
        qe = self.qe
        waves = self.waves
        waves_cm = waves.copy()*10**-7

        # convert the sed to per cm
        n_lambda = sed* c/waves_cm
        # find delta_wave
        dWave = self.find_delta_wavelengths(waves_cm)
        photons = n_lambda*dWave

        photons_mean = scipy.integrate.simps(photons*phi, waves)
        photons = photons_mean * scipy.integrate.simps(qe*waves, waves)
        logN_photons = np.log10(photons)
        return logN_photons

    def find_delta_wavelengths (self, waves) :
        x0 = np.nonzero(waves > 0)[0]
        x1 = x0+1; 
        bad_x = np.nonzero(x1 >= len(waves))[0]
        last_good_x = bad_x[0]-2
        x1[bad_x] = x1[last_good_x]
        delta_waves = waves[x1] - waves[x0]
        delta_waves[bad_x] = delta_waves[last_good_x]
        return delta_waves

    def interpolateSed (self, wavelengths, sed_waves, sed_flux) :
        """Interpolate the SED to the wavelengths giving using a 
        spline that passes through the sed flux points.
        """
        from scipy import interpolate
        spline = interpolate.InterpolatedUnivariateSpline(sed_waves, sed_flux)
        flux = spline(wavelengths)
        return flux

    def loadSystemResponse(self) :
        """ this now loads a dictionary of sys responses into trans
        indexed by ccd number
        """
        if self.camera == "des" : 
            maxCCD = 63
        elif self.camera == "sdss" :
            maxCCD = 2
        else :
            raise Exception("only sdss or des: {}".format(self.camera))
        sysResponseFile = self.sysResponseFile
        data = np.genfromtxt(sysResponseFile, unpack=True)
        trans = dict()
        for ccd in range(1,maxCCD) :
            trans[ccd] = data[ccd]
        waves = data[0]
        return [waves,trans]

    def plot (self, fig, plotSedType) :
        """Flux that the Blanco would see (aperture area in cm = 10^8)

        The SED can be plotted as 
            ergs/hz/sec/cm^2     = "f_nu"
            ergs/A/sec/cm^2      = "f_lambda"
            photons/hz/sec/cm^2  = "photon_nu"
            photons/A/sec/cm^2   = "photon_lambda"

        We work with the object sed, which is always stored as
        photon_lambda. We can convert at will.

        Here we also adopt the notion that what is fundamental isn't
        ergs/sec/cm^2/Hz or ergs/sec/cm^2/A but ergs/sec/cm^2.
        This means that we plot 
            nu F_nu vs log(nu) or lambda F_lambda vs log(lambda)
        For photons,
            nu photons vs log(nu) or lambda photons vs log(lambda)
        Think of this as per fractional (lambda or nu ) interval
        So there are four additional plots:
            ergs/sec/cm^2        = "energy" =  lambda f_lambda vs log(lambda) 
            photons/sec/cm^2     = "photon" = lambda photon_lambda vs log(lambda)
        """
        import matplotlib.ticker
        h = self.h
        c = self.c
        hc = self.hc
        AperCm = 10**8.
        AperCm = 1.
        photon_nu = self.sed
        waves = self.waves
        waves_cm = waves.copy()*10**-7 #; nm->cm
        do_log_wave = False
        if (plotSedType == "photon_nu") :
            flux_label = "photons$_\\nu$ (photons/s/cm$^2$/Hz)"
            flux = photon_nu
        elif (plotSedType == "f_nu") :
            flux_label = "$f_\\nu$ (ergs/s/cm$^2$/Hz)"
            energyFluxHz = photon_nu * hc/waves_cm
            flux = energyFluxHz
        elif (plotSedType == "f_lambda")  :
            flux_label = "$f_\lambda$ (ergs/s/cm$^2$/A)"
            energyFluxCM = photon_nu * h *c**2 / waves_cm**3
            flux = energyFluxCM/AperCm
        elif (plotSedType == "photon_lambda") :
            flux_label = "photons$_\lambda$ (photons/s/cm$^2$/A)"
            photonFluxCM = photon_nu * c / waves_cm**2 
            flux = photonFluxCM/AperCm
        elif (plotSedType == "energy") :
            flux_label = "$\\lambda f_\lambda$ (ergs/s/cm$^2$)"
            energyFlux = waves_cm * photon_nu * h *c**2 / waves_cm**3
            flux = energyFlux/AperCm
            do_log_wave = True
        elif (plotSedType == "photon") :
            flux_label = "$\\lambda$ photons$_\lambda$ (photons/s/cm$^2$)"
            photonFlux = waves_cm * photon_nu * c / waves_cm**2 
            flux = photonFlux/AperCm
            do_log_wave = True
        else :
            raise NameError('sed type not known: {}'.format(plotSedType))

        qe = self.qe
        waves = self.waves

        max = flux.max()
        max_q = qe.max()
        qe = 1.0*(max/max_q)* qe

        ax = fig.add_subplot(1,1,1)
        plt.plot(waves,flux,color="b")
        plt.plot(waves,qe,color="k")
        if do_log_wave :
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(300,1100)
        plt.xlabel("wavelength (nm)")
        plt.ylabel(flux_label)


