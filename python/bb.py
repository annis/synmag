import synMagnitude 
import numpy as np

            # this one keeps interesting NS scaling
            # L in units of 10^{40} ergs/sec, d in units of 100 Mpc, T in units of 1000 K
            # so the above scaling is
            #T=2.255; d= 1.0; L = 1.38
            #flux = flux * 4.4e-22  * L/d**2/T**4
def setup(temperature, fluxScaling = 4.4e-22, L40 = 1, distance = 100.0) :
    t1000 = temperature/1000.
    d100 = distance/100.
    fluxScaling =fluxScaling * L40 / d100**2 / t1000**4
    text = "/home/s1/annis/daedalean/synmag/sysresp/DES-{}.txt"
    decam,mag=dict(),dict()
    for f in ["g","r","i","z","y"] :
        decam[f] = synMagnitude.synMag(text.format(f)); 
        decam[f].set_ccd(36)
        decam[f].set_blackbody_sed(temperature, fluxScaling=fluxScaling)
        mag[f] = decam[f].magnitude()
    return decam, mag

def energy (decam, temperature, filter) :
    nmToCm = 1e-7
    h = 6.63e-27 ;# erg sec
    c = 3.00e10 ;# cm/s
    sigma_b = 5.67e-5 ;# stefan-boltzman constant,  erg /cm^2/sec/K^4

    # the 1/pi is due to the Int B_nu dnu = sigma T^4/pi
    bb_total_energy_flux_unit_area = sigma_b*temperature**4/np.pi
    print "log(bb_total_energy_flux_unit_area) ", np.log10(bb_total_energy_flux_unit_area) , " ergs/s/cm^2"

    Lambda_e = decam[filter].Lambda_e()
    Energy_e = h*c / (Lambda_e * nmToCm)
    print "effective wavelength: ", Lambda_e, "nm"
    print "log(effective energy per photon): ",np.log10(Energy_e), " ergs"

    log_photon_flux = decam[filter].logN_photons_mean(); 
    print "log(effective number of photons) :", log_photon_flux, " photons/s/cm^2"
    
    energy = Energy_e * 10**log_photon_flux
    print "filter total energy: ", energy, " ergs/s/cm^2"

    fraction_of_blackbody = energy/bb_total_energy_flux_unit_area; 
    print "fraction of blackbody energy: ", fraction_of_blackbody 
