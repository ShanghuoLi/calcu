import io
import os
import numpy as np

# import local things
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
        print ("frequency = {} GHz".format(self.freq))
        # dust opacity (Hildebrand 1983): K_v =10*(V/1.2THz)^beta  here we using beta=1.5
        K_nu = 10*(self.freq/1200.)**self.beta
        print ("K_nu= {} cm^2 g^-1".format(K_nu) ) #cm^2/g)
        #Planck function; 2h nu**3/c**2 * 1.0 / exp(h nu/ kT) -1.0
        b1=2.0*h*( (self.freq*1e9)**3.0)/(c**2.0)

        Pixelsize = 3600. * abs(self.header['CDELT1'])
        PixelAreaArcsec = 3600. * abs(self.header['CDELT1']) * 3600. * abs(self.header['CDELT2'])
        print ('Pixel size:{} arcsec, PixelAreaArcsec={} arcsec^2'.format(Pixelsize,PixelAreaArcsec))
        bmaj = self.header['BMAJ'] * 3600
        bmin = self.header['BMIN'] * 3600
        print ("beam size: bmaj x bmin = {} arcsec x {} arcsec".format( bmaj, bmin ))
        sigma2fwhm = (8*np.log(2))**0.5
        beamArea = 2*np.pi*bmaj*bmin/(sigma2fwhm**2)
        numpixbeam = beamArea/PixelAreaArcsec
        print ("number of pixel in a beam = {}".format(numpixbeam) )
        effectivebeam = ( beamArea/(2*np.pi/(sigma2fwhm**2)) )**0.5
        print ("effectivebeam(FWHM) ={} in arcsec ".format(effectivebeam))
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
    flux: a float/int value, numpy.ndarray or 2D fits image
        A 2D array of the data to be analyzed.
        If the unit is Jy/beam*km/s, the Jyperbeam2K need to be provided.
        If the unit is K*km/s, the Jyperbeam2K=1 (default).  

    name: str
        line name, e.g., CO, 13CO, H13CO+ ....

    J: str
        quantum numbers, e.g., '1-0', '2-1', '3-2'. 

    Tex: numpy.ndarray, a float or int value
        excitation temperature (K).

    Jyperbeam2K: float
        Jy/beam to K;
        if the unit of the flux is Jy/beam*km/s, the convertion factor need to be provided.
        if the unit of flux is K*km/s, the Jyperbeam2K=1 (default)

    tau: float or int
        optical depth, default is 0, no optical depth correction. 


    -------------------- output --------------------
    column: numpy.ndarray or a value
        column density (cm^-2). 

    Sijmu2: Sij*mu^2 (D^2). 

    B0: rigrid rotor rotattion constant (Hz). 

    freq: freqeency (Hz). 

    Eu: upper energy (K). 

    '''
    ##---- define functions 

    def __init__(self, flux=None, pb=None, name='CO', J='1-0', Tex=15., Jyperbeam2K=1, tau=0):
        self.name = name
        self.Tex = Tex
        self.Jyperbeam2K = Jyperbeam2K
        self.tau = tau
        
        ##----- load the flux ------
        if isinstance(flux, (int, float)):
            Iflux = flux
        elif (isinstance(flux, np.ndarray)) and (flux.ndim==2):
            Iflux = flux
        elif flux[-5:]=='.fits':
            Inputflux = fits.open(flux)[0]
            if Inputflux.header['NAXIS'] == 2:
                Iflux = Inputflux.data
            elif Inputflux.header['NAXIS'] > 2:
                raise TypeError("flux should be given as 2D fits image.")
        else:
            raise TypeError("flux should be given as a single int or float value, 2D numpy.array, or 2D fits.")
        


        ##----- load the primary beam correction factor ------
        if pb is None:
            pbcor = 1
            # print("No primary beam correction is applied to the datacube")
        elif isinstance(pb, (int, float)):
            pbcor = pb
        elif (isinstance(pb, np.ndarray)) and (pb.ndim==2):
            pbcor = pb
        elif pb[-5:]=='.fits':
            Inputpb = fits.open(pb)[0]
            if Inputpb.header['NAXIS'] == 2:
                pbcor = Inputpb.data
            elif Inputpb.header['NAXIS'] == 4:
                pbdata = Inputpb.data[0]
                for i in range( Inputpb.header['NAXIS1'] ):
                    if np.isnan( np.nanmax(pb[i,:,:]) ) == False:
                        pbcor = pbdata[i,:,:]
                        break            
            elif Inputpb.header['NAXIS'] == 3:
                pbdata = Inputpb.data
                for i in range( Inputpb.header['NAXIS1'] ):
                    if np.isnan( np.nanmax(pbdata[i,:,:]) ) == False:
                        pbcor = pbdata[i,:,:]
                        break
            elif InputMap.header['NAXIS'] > 4:
                raise TypeError("pb should be given as 2D image or 3D cube.")
            else:
                raise TypeError("pb should be given as 2D image or 3D cube.")
        else:
            raise TypeError("pb should be given as a single int or float value, 2D numpy.array, or 2D or 3D fits.")

        self.pb = pbcor
        if isinstance(pbcor, np.ndarray):
            if pbcor.shape != Iflux.shape:
                raise TypeError("'pb' shape is not match to flux image")

        self.flux = (Iflux/pbcor) * self.Jyperbeam2K


        if name=='CO':
            self.molecule = line_CO
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs )) 
            self.note = self.molecule['note']
            # self.mu = mu_CO           
            # Qns   = self.molecule['QNs'][J]
            # gu    = self.molecule['gu'][J]
            freq   = self.molecule['freq'][J] * 1e9  # GHz to HZ
            Eu     = self.molecule['Eu'][J]     # K
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36  # D^2 to esu^2 cm^2
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6  # MHz to Hz

            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (k*self.Tex/(h*B0) + 1./3)


        ##------- H13COP 
        elif name=='H13CO+':
            # self.mu = mu_H3COp
            self.molecule = line_H13COp
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs )) 
            self.note = self.molecule['note']
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (k*self.Tex/(h*B0) + 1./3)


        ##------- 13CO 
        elif name=='13CO':
            self.molecule = line_13CO
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs )) 
            self.note = self.molecule['note']
            # self.mu = mu_CO           
            # Ju    = self.molecule['QNs'][J]
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
            self.molecule = line_C17O
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs )) 
            self.note = self.molecule['note']
            # self.mu = mu_CO           
            # Ju    = self.molecule['QNs'][J]
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
            self.molecule = line_C18O
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs )) 
            self.note = self.molecule['note']
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
            self.molecule = line_SiO
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            # Ju    = self.molecule['QNs'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (k*self.Tex/(h*B0) + 1./3)

            
        ##------- CS
        elif name=='CS':
            self.molecule = line_CS
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            # Ju    = self.molecule['QNs'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (1.85*self.Tex - 3.32)


        ##------- C34S
        elif name=='C34S':
            self.molecule = line_C34S
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            # Ju    = self.molecule['QNs'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (1.85*self.Tex - 3.32)


        ##------- 13CS
        elif name=='13CS':
            self.molecule = line_13CS
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            # Ju    = self.molecule['QNs'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (1.85*self.Tex - 3.32)



        ##------- N2H+
        elif name=='N2H+':
            self.molecule = line_N2Hp
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            # Ju    = self.molecule['QNs'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (0.9*(self.Tex)**(1.28) + 37.18)


######-------------------------------------------------------------------------------------
### duterium-molecules
######-------------------------------------------------------------------------------------
        ##------- NH2D ----------------------
        elif name=='NH2D':
            # self.mu = mu_NH2D
            # B0 = B0_NH2D
            # Ju = line_NH2D['QNs'][J]
            self.molecule = line_NH2D
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            # gu     = self.molecule['gu'][J]
            freq   = self.molecule['freq'][J] * 1e9 
            Eu     = self.molecule['Eu'][J]
            Ri     = self.molecule['Ri'][J]
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (0.751*(self.Tex)**(3/2) + 3.899)


        ##------- N2D+
        elif name=='N2D+':
            self.molecule = line_N2Dp
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            # Ju    = self.molecule['QNs'][J]
            freq   = self.molecule['freq'][J] * 1e9
            # gu    = self.molecule['gu'][J]
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (4.87*self.Tex + 2.81)



        ##------- DCN
        elif name=='DCN':
            self.molecule = line_DCN
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (0.1*(self.Tex)**(3/2) + 50.51)


        ##------- DCO+
        elif name=='DCO+':
            self.molecule = line_DCOp
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (0.58*self.Tex + 0.34)


        ##------- C2D
        elif name=='C2D':
            self.molecule = line_C2D
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (3.47*self.Tex + 2.06)



        ###-------------------------------------------------------------------------------
        ##  slightly-asymmetric rotor molecule
        ###-------------------------------------------------------------------------------

        #------- H2CO
        elif name=='H2CO':
            self.molecule   = line_H2CO
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            # self.note = self.molecule['note']
            freq  = self.molecule['freq'][J] * 1e9
            Eu    = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri    = self.molecule['Ri'][J]
            A0    = self.molecule['B0'][0] * 1e6
            B0    = self.molecule['B0'][1] * 1e6
            C0    = self.molecule['B0'][2] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (1./2) * ( (np.pi*(k*self.Tex)**3 ) / (h**3 * A0*B0*C0) )**0.5

        #------- CH3OH
        elif name=='CH3OH':
            self.molecule   = line_H2CO
            if J in self.molecule['QNs']:
                self.QNs = J
                J = self.molecule['QNs'].index(J)
            else:
                raise TypeError("please provide a correct Quantum numbers (transition).")
            print("---->>> Line {} {}  <<<-------".format(name, self.QNs ))
            self.note = self.molecule['note']
            freq  = self.molecule['freq'][J] * 1e9
            Eu    = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri    = self.molecule['Ri'][J]
            A0    = self.molecule['B0'][0] * 1e6
            B0    = self.molecule['B0'][1] * 1e6
            C0    = self.molecule['B0'][2] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (1./3) * ( (np.pi*(k*self.Tex)**3 ) / (h**3 * A0*B0*C0) )**0.5




        ###-------------------------------------------------------------------------------
        ###-------------------------------------------------------------------------------

        Ncolum = const0 * Qrot * (np.exp(Eu/self.Tex) / (np.exp(h*freq/(k*self.Tex))-1)) * \
        (1/(J_nu(self.Tex,freq) - J_nu(Tbg,freq))) * self.flux *1e5 

        ### optical depth correction
        if isinstance(tau, (int, float)):
            self.tau = tau
            if self.tau==0:
                self.column = Ncolum
            else:
                self.column = Ncolum * ( self.tau/(1 - np.exp(-self.tau) ) )
        elif (isinstance(flux, np.ndarray)) and (flux.ndim==2):
            self.tau = tau
            self.column = Ncolum * ( self.tau/(1 - np.exp(-self.tau) ) )
        elif tau[-5:]=='.fits':
            Inputtau = fits.open(tau)[0]
            if Inputtau.header['NAXIS'] == 2:
                self.tau = Inputtau.data
                self.column = Ncolum * ( self.tau/(1 - np.exp(-self.tau) ) )
            elif Inputtau.header['NAXIS'] > 2:
                raise TypeError("tau should be given as 2D fits image.")
        else:
            raise TypeError("tau (optical depth) should be given as a single int or float value, 2D numpy.array, or 2D fits.")
        ###-------------------------------------------------------------------------------
        ###-------------------------------------------------------------------------------       
        # 1e5 is the factor of km/s to cm/s
        # const0 = 3 * k / (8 *(np.pi)**3 * Sijmu2 * Ri *freq)
        # Qrot = (k*self.Tex/(h*B0) + 1./3)
        # self.column = const0 * (Qrot) * \
        #     (np.exp(Eu/self.Tex) / (np.exp(h*freq/(k*self.Tex))-1)) * \
        #     (self.Tex/(self.Tex - Tbg)) * \
        #     self.flux *1e5  

        ###-------------------------------------------------------------------------------
        ###-------------------------------------------------------------------------------
