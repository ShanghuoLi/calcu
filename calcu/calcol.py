import io
import os
import numpy as np
from .line_info import *


###------------------------------------------------------------------------
###------------------------------------------------------------------------
##
## To calculate the H2 column density and molecular column density
##
###------------------------------------------------------------------------
###------------------------------------------------------------------------
##---------------
## Constants
##____________________________
h = 6.634e-27        #erg s; Planck's constant
k = 1.38e-16         # erg/K; Boltzmann's constant
c = 3e10             #cm/s, light speed
Msolar = 1.98892e33  #gram; solar mass
mH = 1.673534e-24    #gram; hydrogen atom mass
kpc2cm = 3.0857e21   #kpc to cm
ug = 2.8            # mean molecular weight of the ISM
co2h2 = 1e4         # H2 to CO ratio

Tbg = 2.73          ## K, background temperature; CMB
###------------------------------------------------------------------------


##------------------------------------------------------------------------
##------------ using the dust continuum to calculate the H2 column density
class calcol_cont(object):
    '''
    using the dust continuum to calculate the H2 column density

    Parameters 
    -------------- input ----------
    header: astropy.io.fits.Header 
        The header corresponding to the image array

    flux: numpy.ndarray or a value.
        A 2d array of the data (unit is Jy/beam).

    Td: dust temperature (K)

    eta: gas-to-dust ratio 

    beta: dust beta 

    -------------- output ----------
    column: column density (cm^-2)

    freq: frequency (GHz)

    effbeam: effective beam (arcsec)

    '''
    def __init__(self, header, flux=1, Td=15, eta=100, beta=1.5):
        self.flux = flux # unit is Jy/beam
        self.Td = Td     # dust temperature
        self.beta = beta    #
        self.eta = eta      # the gas-to-dust ratio
        self.header = header

        self.freq = self.header['RESTFRQ'] * 1e-9  # GHz --> Hz
        print ("frequency = %s GHz"%(self.freq))
        # dust opacity (Hildebrand 1983): K_v =10*(V/1.2THz)^beta  here we using beta=1.5
        K_nu = 10*(self.freq/1200.)**self.beta
        print ("K_nu= %s cm^2 g^-1"%(K_nu) ) #cm^2/g)
        #Planck function; 2h nu**3/c**2 * 1.0 / exp(h nu/ kT) -1.0
        b1=2.0*h*( (self.freq*1e9)**3.0)/(c**2.0)

        Pixelsize = 3600. * abs(self.header['CDELT1'])
        PixelAreaArcsec = 3600. * abs(self.header['CDELT1']) * 3600. * abs(self.header['CDELT2'])
        print ('Pixel size:%s arcsec, PixelAreaArcsec=%s arcsec^2'%(Pixelsize,PixelAreaArcsec))
        bmaj = self.header['BMAJ'] * 3600
        bmin = self.header['BMIN'] * 3600
        print ("beam size: bmaj x bmin = %s arcsec x %s arcsec"%( bmaj, bmin ))
        sigma2fwhm = (8*np.log(2))**0.5
        beamArea = 2*np.pi*bmaj*bmin/(sigma2fwhm**2)
        numpixbeam = beamArea/PixelAreaArcsec
        print ("number of pixel in a beam = %s"%(numpixbeam) )
        effectivebeam = ( beamArea/(2*np.pi/(sigma2fwhm**2)) )**0.5
        print ("effectivebeam(FWHM) =%s in arcsec "%(effectivebeam))
        self.effbeam = effectivebeam
        
        # Using the dust emission to calculate the H2 column density
        self.column = 5.43635*1e29 * self.eta * self.flux * (self.freq**-3) * \
        (self.effbeam**-2) * (np.exp( 0.048 * (self.freq) * (self.Td)**-1 ) - 1) * \
        ((self.freq/1200.)**-1.5)
#         err_N_H2_cont = 5.43635*10**29 * self.eta * rms * (self.freq**-3) * (effectivebeam**-2) *\
#         (np.exp( 0.048 * (self.freq) * (self.Td)**-1 ) - 1) * ((self.freq/1200.)**-1.5)
###------------------------------------------------------------------------
###------------------------------------------------------------------------



##--------------------------------------------------------------------------------------
## define function
##--------------------------------------------------------------------------------------

def J_nu(T, freq):
    J_nu = (h*freq/k) / ( np.exp(h*freq/(k*T)) - 1 )
    return J_nu

# def Nfun(const0, Qrot, Eu,  Tex, freq, flux):
#     Nfun = const0 * Qrot * (np.exp(Eu/Tex) / (np.exp(h*freq/(k*Tex))-1)) * \
#     (1/(J_nu(Tex,freq) - J_nu(Tbg,freq))) * flux *1e5  
#     return Nfun
# 1e5 is the factor of km/s to cm/s

##---------------------
# class column(object):
class calcol_line(object):
    '''
    using the line emission to calculate the total molecular column density

    Parameters
    -------------------- input -------------------- 
    flux: numpy.ndarray or a value
        A 2D array of the data to be analyzed.
        If the unit is Jy/beam*km/s, the Jyperbeam2K need to be provided.
        If the unit is K*km/s, the Jyperbeam2K=1 (default).  

    name: line name
        e.g., CO, 13CO, H13CO+ ....

    J: rotational quantum number of the upper state, 
        e.g., J=1 mean transition "1-0". 

    Tex: numpy.ndarray or a value
        excitation temperature (K).

    Jyperbeam2K: Jy/beam to K
        if the unit of the flux is Jy/beam*km/s, the convertion factor need to be provided.
        if the unit of flux is K*km/s, the Jyperbeam2K=1 (default)

    tau: 
    optical depth 
        default is 0, no optical depth correction. 


    -------------------- output --------------------
    column: column density (cm^-2). 

    Sijmu2: Sij*mu^2 (D). 

    B0: rigrid rotor rotattion constant (Hz). 

    freq: freqeency (Hz). 

    Eu: upper energy (K). 

    '''
    ##---- define functions 

    def __init__(self, flux=1, name='CO', J=1, Tex=15, Jyperbeam2K=1, tau=0):
        self.name = name
        self.J = J
        J = J-1
        self.Tex = Tex
        self.Jyperbeam2K = Jyperbeam2K
        self.flux = flux * self.Jyperbeam2K
        self.tau = tau

        if name=='CO':
            print("---->>> Line 12CO %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_CO
            # print(self.molecule['note'])
            # self.mu = mu_CO           
            # Ju    = self.molecule['Ju'][J]
            freq   = self.molecule['freq'][J] * 1e9  # GHz to HZ
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]     # K
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36  # D^2 to esu^2 cm^2
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6  # MHz to Hz

            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (k*self.Tex/(h*B0) + 1./3)


        ##------- H13COP 
        elif name=='H13CO+':
            print("---->>> Line H13CO+ %s-%s  <<<-------"%(J+1, J))
            # self.mu = mu_H3COp
            self.molecule = line_H13COp
            # print(self.molecule['note'])
            # Ju    = self.molecule['Ju'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (k*self.Tex/(h*B0) + 1./3)


        ##------- 13CO 
        elif name=='13CO':
            print("---->>> Line 13CO %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_13CO
            # print(self.molecule['note'])
            # self.mu = mu_CO           
            # Ju    = self.molecule['Ju'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6

            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (k*self.Tex/(h*B0) + 1./3)


        ##------- C17O 
        elif name=='C17O':
            print("---->>> Line C17O %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_C17O
            # print(self.molecule['note'])
            # self.mu = mu_CO           
            # Ju    = self.molecule['Ju'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6

            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (2.23*self.Tex + 1.92)


        ##------- C18O 
        elif name=='C18O':
            print("---->>> Line C18O %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_C18O
            # print(self.molecule['note'])
            # self.mu = mu_CO           
            # Ju    = self.molecule['Ju'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6

            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (k*self.Tex/(h*B0) + 1./3)


 
        ##------- SiO 
        elif name=='SiO':
            print("---->>> Line SiO %s-%s  <<<-------"%(J+1, J))
            # self.mu = mu_SiO 
            self.molecule = line_SiO
            # print(self.molecule['note'])
            # Ju    = self.molecule['Ju'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (k*self.Tex/(h*B0) + 1./3)


        ###-------------------------------------------------------------------------------
        ###-------------------------------------------------------------------------------

        Ncolum = const0 * Qrot * (np.exp(Eu/self.Tex) / (np.exp(h*freq/(k*self.Tex))-1)) * \
        (1/(J_nu(self.Tex,freq) - J_nu(Tbg,freq))) * self.flux *1e5 

        if self.tau==0:
            self.column = Ncolum
        else:
            self.column = Ncolum * ( self.tau/(1 - np.exp(-self.tau) ) )
        # 1e5 is the factor of km/s to cm/s
        # const0 = 3 * k / (8 *(np.pi)**3 * Sijmu2 * Ri *freq)
        # Qrot = (k*self.Tex/(h*B0) + 1./3)
        # self.column = const0 * (Qrot) * \
        #     (np.exp(Eu/self.Tex) / (np.exp(h*freq/(k*self.Tex))-1)) * \
        #     (self.Tex/(self.Tex - Tbg)) * \
        #     self.flux *1e5  

        ###-------------------------------------------------------------------------------
        ###-------------------------------------------------------------------------------
