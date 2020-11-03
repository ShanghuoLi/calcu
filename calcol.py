import io
import os
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from line_info import *

##--------------------------------------------------------------------------------------
## define function
##--------------------------------------------------------------------------------------
class immoment0(object):
    ''' 
    Calculate the moment-0
    Parameters
    -------------- input ---------- 
    fname: file name            
    vl: a number-like object
        staring velocity in km/s     
    vu: a number-like object
        ending velocity in km/s   
    pb: numpy.ndarray (2D or 3D) or a number-like object
        primary beam correction;
    -------------- output ---------- 
    momt0: numpy.ndarray
        velocity-integrated intensity map
    '''
    def __init__(self, fname, vl=0, vu=1, pb=None):
        self.vl = vl
        self.vu = vu
        if pb[-5:] =='.fits':
            InputMap = fits.open(fname)[0]
        else:
            raise TypeError("data should be given as 3D-spectrum cube fits file.")

        if InputMap.header['NAXIS'] == 4:
            img = InputMap.data[0]
        elif InputMap.header['NAXIS'] == 3:
            img = InputMap.data
        elif (InputMap.header['NAXIS'] < 3) or (InputMap.header['NAXIS'] > 4):
            raise TypeError("data should be given as 3D-spectrum cube fits file.")
        else:
            raise TypeError("data should be given as 3D-spectrum cube fits file.")
        self.onevpix = np.abs(InputMap.header['CDELT3']) * 0.001  # from m/s to km/s
        v0 = InputMap.header['CRVAL3'] * 0.001  # from m/s to km/s, reference velocity
        v0pix = InputMap.header['CRPIX3']       # channel width in velocity

        self.vl_chan = np.int( np.abs((self.vl - v0)) / np.abs(self.onevpix) + v0pix )
        self.vu_chan = np.int( np.abs((self.vu - v0)) / np.abs(self.onevpix) + v0pix )
        print ( "----------------------" )
        print ( 'starting velocity', self.vl,'km/s at channel', self.vl_chan )
        print ( 'ending velocity  ', self.vu,'km/s at channel', self.vu_chan )

        if self.vl_chan < self.vu_chan:
            vch_start,vch_stop = self.vl_chan,self.vu_chan
        elif self.vl_chan > self.vu_chan:
            vch_start,vch_stop = self.vu_chan,self.vl_chan

        momt0 = np.zeros((img.shape[1],img.shape[2]))
        momt0 = np.sum(img[vch_start:vch_stop+1,:,:], axis=0) * np.abs(onevpix) # Sdv# JY/beam * km/s
        
        if pb==None:
            self.momt0 = momt0
        elif isinstance(pb, (int, float)):
            self.momt0 = momt0/pb
        elif isinstance(pb, np.ndarray) and (len(pb.shape)==2):
            self.momt0 = momt0/pb
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
            self.momt0 = momt0/pbcor
        else:
            raise TypeError("pb should be given as a single int or float value, 2D numpy.array, or 2D or 3D fits.")
##--------------------------------------------------------------------------------------




##--------------------------------------------------------------------------------------
class headerinfo(object):
    '''
    Get some useful information from header

    input is the file name 
    -------------- output ----------
    header: the header of the file   

    onevpix: the velocity speration of a channel (deltav) 

    v0: the first channel velocity 

    v0pix: Pixel coordinate of reference point for velocity

    onepix: the pixel size 

    PixelAreaArcsec: pixel area in arcsec^2 

    bmaj: the beam major axis 

    bmin: the beam minor axis 

    freq: the reference freqeuncy 

    Jyperbeam2K: a convertion factor of Jy/beam to K  

    hduwcs2d: the header corresponing to the 2D image
    '''
    def __init__(self, fname):
        InputMap = fits.open(fname)[0]
        self.header = InputMap.header
        self.onevpix = InputMap.header['CDELT3'] * 0.001  # from m/s to km/s
        self.v0 = InputMap.header['CRVAL3'] * 0.001  # from m/s to km/s
        self.v0pix = InputMap.header['CRPIX3']
        self.onepix = 3600. * abs(InputMap.header['CDELT1'])  #pixel size in arcsec
    
        self.PixelAreaArcsec = 3600. * abs(InputMap.header['CDELT1']) * \
                            3600. * abs(InputMap.header['CDELT2'])
        self.bmaj = InputMap.header['BMAJ'] * 3600
        self.bmin = InputMap.header['BMIN'] * 3600
        print ( "major_beam=%s minor_beam=%s"%(self.bmaj, self.bmin) )
        self.freq = InputMap.header['RESTFRQ'] / 1e9     # GHz

        if InputMap.header['BUNIT'] == 'beam-1 Jy':
            self.Jyperbeam2K = 1.222e6/(self.freq**2 * self.bmaj * self.bmin)  # freq in GHz, bmaj and bmin in arcsec
            print ( "Jyperbeam_to_K = %.3f"%(self.Jyperbeam2K) )
        elif InputMap.header['BUNIT'] == 'K':
            self.Jyperbeam2K = 1
            print("The unit of datacube is K")
            print ( "Jyperbeam_to_K = %.3f"%(self.Jyperbeam2K) )

        hduhead = self.header
        if ('PC01_01' in hduhead) and ('PC01_04' in hduhead):
            rm_key = ['NAXIS4','CRPIX4','CDELT4', 'CUNIT4', 'CTYPE4', 'CRVAL4',
            'NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
            'PC01_01','PC02_01','PC03_01','PC04_01',
            'PC01_02','PC02_02','PC03_02','PC04_02',
            'PC01_03','PC02_03','PC03_03','PC04_03',
            'PC01_04','PC02_04','PC03_04','PC04_04']
        elif ('PC1_1' in hduhead) and ('PC1_4' in hduhead):
            rm_key = ['NAXIS4','CRPIX4','CDELT4', 'CUNIT4', 'CTYPE4', 'CRVAL4',
            'NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
            'PC1_1', 'PC1_2', 'PC1_3', 'PC1_4', 
            'PC2_1', 'PC2_2', 'PC2_3', 'PC2_4',
            'PC3_1', 'PC3_2', 'PC3_3', 'PC3_4',
            'PC4_1', 'PC4_2', 'PC4_3', 'PC4_4']

        elif ('PC01_01' in hduhead) and ('PC01_03' in hduhead):
            rm_key = ['NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
            'PC01_01','PC02_01','PC03_01',
            'PC01_02','PC02_02','PC03_02',
            'PC01_03','PC02_03','PC03_03']

        elif ('PC1_1' in hduhead) and ('PC1_3' in hduhead):
            rm_key = ['NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
            'PC1_1', 'PC1_2', 'PC1_3', 
            'PC2_1', 'PC2_2', 'PC2_3',
            'PC3_1', 'PC3_2', 'PC3_3']
        elif ('CRPIX3' in hduhead):
            rm_key = ['NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
        else:
            rm_key = None

        if rm_key != None:
            for key_i in rm_key:
                hduhead.remove(key_i)
        hduhead['NAXIS'] = 2
        hduhead['WCSAXES'] = 2
        self.hduhead = hduhead
        self.hduwcs2d = wcs.WCS(hduhead)

##--------------------------------------------------------------------------------------


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


#------------------------------------------------------
#------------------------------------------------------
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



        ##------- 13CS
        elif name=='13CS':
            print("---->>> Line 13CS %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_13CS
            # print(self.molecule['note'])
            # Ju    = self.molecule['Ju'][J]
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
            print("---->>> Line N2D+ %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_N2Hp
            # print(self.molecule['note'])
            # Ju    = self.molecule['Ju'][J]
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
            print("---->>> Line NH2D (1(1,1)-1(0,1))  <<<-------")
            # self.mu = mu_NH2D
            # B0 = B0_NH2D
            # Ju = line_NH2D['Ju'][J]
            self.molecule = line_NH2D
            # print(self.molecule['note'])
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            # gu     = self.molecule['gu'][J]
            freq   = self.molecule['freq'][J] * 1e9 
            Eu     = self.molecule['Eu'][J]
            Ri     = self.molecule['Ri'][J]
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (0.751*(self.Tex)**(3/2) + 3.899)


        ##------- N2D+
        elif name=='N2D+':
            print("---->>> Line N2D+ %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_N2Dp
            # print(self.molecule['note'])
            # Ju    = self.molecule['Ju'][J]
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
            print("---->>> Line DCN %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_DCN
            # print(self.molecule['note'])
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot = (0.1*(self.Tex)**(3/2) + 50.51)


        ##------- DCO+
        elif name=='DCO+':
            print("---->>> Line DCO+ %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_DCOp
            # print(self.molecule['note'])
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (0.58*self.Tex + 0.34)


        ##------- C2D
        elif name=='C2D':
            print("---->>> Line C2D %s-%s  <<<-------"%(J+1, J))
            self.molecule = line_C2D
            # print(self.molecule['note'])
            freq   = self.molecule['freq'][J] * 1e9
            Eu     = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri     = self.molecule['Ri'][J]
            B0     = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (3.47*self.Tex + 2.06)





        # ##------- C2D
        # elif name=='C2D':
        #     print("---->>> Line C2D %s-%s  <<<-------"%(J+1, J))
        #     self.molecule = line_C2D
        #     pri          #  'B0':[281970.56,38833.987,34004.244],
            # mu2 = self.molecule['Sijmu2'][J] * 1e-36
        #     Ri    = self.molecule['Ri'][J]
        #     B0    = self.molecule['B0'][0] * 1e6
        
        #     const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
        #     Qrot = (3.47*self.Tex + 2.06)





        #------- H2CO
        elif name=='H2CO':
            self.molecule   = line_H2CO
            print("---->>> Line H2CO %s  <<<-------"%(self.molecule['Ju'][J]))
            # # print(self.molecule['note'])
            freq  = self.molecule['freq'][J] * 1e9
            Eu    = self.molecule['Eu'][J]
            Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
            Ri    = self.molecule['Ri'][J]
            # B0    = self.molecule['B0'][0] * 1e6
        
            const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
            Qrot   = (0.56*(self.Tex)**(3./2) - 0.64)


        # #------- CH3OH
        # elif name=='CH3OH':
        #     self.molecule   = line_H2CO
        #     print("---->>> Line CH3OH %s  <<<-------"%(self.molecule['Ju'][J]))
        #     # print(self.molecule['note'])
        #     freq  = self.molecule['freq'][J] * 1e9
        #     Eu    = self.molecule['Eu'][J]
        #     Sijmu2 = self.molecule['Sijmu2'][J] * 1e-36
        #     Ri    = self.molecule['Ri'][J]
        #     # B0    = self.molecule['B0'][0] * 1e6
        
        #     const0 = 3 * h / (8 *(np.pi)**3 * Sijmu2 * Ri )
        #     Qrot   = (4.58*(self.Tex)**(3./2) - 70.26)




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
