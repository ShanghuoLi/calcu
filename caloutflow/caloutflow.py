from __future__ import division
import io
import os
import warnings
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from calcol import calcol_line

##----------------------------
## Constants
##____________________________
h = 6.634e-27        #erg s
k = 1.38e-16         # erg/K
c = 3e10             #cm/s, light speed
Msolar = 1.98892e33  # gram
mH = 1.673534e-24    #gram
kpc2cm = 3.0857e21   #kpc to cm
ug = 2.8             # mean molecular weight of the ISM; See Kauffmann et al. 2008, A&A 487, 993 (http://dx.doi.org/10.1051/0004-6361:200809481)
co2h2 = 1e-4         # CO to H2 ratio
Tbg = 2.73           ## K, background temperature; CMB
###------------------------------------------------------------------------


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

        if fname[-5:] =='.fits':
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
        self.onevpix = np.abs(InputMap.header['CDELT3']) * 0.001  # from m/s to km/s; channel width in velocity
        v0 = InputMap.header['CRVAL3'] * 0.001  # from m/s to km/s, the first channel velocity  velocity
        v0pix = InputMap.header['CRPIX3']       #Pixel coordinate of reference point for velocity

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
        momt0 = np.sum(img[vch_start:vch_stop+1,:,:], axis=0) * np.abs(self.onevpix) # Sdv# JY/beam * km/s

        if pb is None:
            self.momt0 = momt0
        elif isinstance(pb, (int, float)):
            self.momt0 = momt0/pb
        elif (isinstance(pb, np.ndarray)) and (pb.ndim==2):
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


##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
def h2mass(N_mol=None, D=None, onepix=None, Xmol=None, decl=0.):
    ##--------
    if (D is not None) and (onepix is not None):
        xx = D * kpc2cm * onepix * (np.pi/(180. * 60 * 60 )) * (np.cos(np.pi * (abs(decl)/180.))) #cm
        yy = D * kpc2cm * onepix * (np.pi/(180. * 60 * 60 ))   #cm
    else:
        raise TypeError("'distance' should be provided to calculate the solid angle.")
    ##--------
    if (N_mol is not None) and (Xmol is not None): 
        N_H2 = N_mol / Xmol        ## molecular column density to H2 column density
    elif N_mol is None:
        raise TypeError("The total column density of molecule should be provided")
    elif Xmol is None:
        raise TypeError("The fractional abundance of molecule relative to H2 should be provided.")
    else:
        raise TypeError("The molecular column density and its fractional abundance should be provided to calculate the H2 column density.")
    N2_mass = N_H2 * ug* mH * xx * yy   # gram
    N2_mass = N2_mass / Msolar         # Msolar
    return N2_mass
##--------------------------------------------------------------------------------------


##--------------------------------------------------------------------------------------
class calf(object):
    '''
    input is the file name 
    -------------- output ----------
    header: the header of the file   

    onevpix: the velocity speration (deltav) 

    v0: the channel velocity 

    onepix: the pixel size 

    PixelAreaArcsec: pixel area in arcsec^2 

    bmaj: the beam major axis 

    bmin: the beam minor axis 

    freq: the reference freqeuncy 

    Jyperbeam2K: a convertion factor of Jy/beam to K  
    '''

    def __init__(self, fname=None, mask=None, pb=None, D=None, Tex=None, vlsr=None, vl=None, vu=None, line=None, J=None, Xmol=None, decl=0., tau=0.):
        self.decl = decl
        self.tau = tau
        self.Xmol = Xmol
        
        if line is not None:
            self.line = line
        else:
            raise TypeError("the name of molecule should be provided.")
        if J is not None:
            self.J = J
        else:
            raise TypeError("rotational quantum number of the upper state should be provided.")

        if D is not None:
            self.D = D
        else:
            raise TypeError("'distance' should be provided to calculate the solid angle.")

        if vlsr is None:
            raise TypeError("vlsr (systematic velocity) should be provided to calculate the solid angle.")
        
        if vl is None:
            raise TypeError("vl (starting velocity) should be provided to calculate the solid angle.")

        if vu is None:
            raise TypeError("vu (ending velocity) should be provided to calculate the solid angle.")

        if Tex is None:
            raise TypeError("Tex should be provided.")
        elif isinstance(Tex, (int, float)):
            self.Tex = Tex
        elif (isinstance(Tex, np.ndarray)) and (Tex.ndim == 2):
            self.Tex = Tex
        elif Tex[-5:]=='.fits':
            InputTex = fits.open(Tex)[0]
            if InputTex.data.ndim == 2:
                self.Tex = InputTex.data
            else:
                raise TypeError("Tex should be given as 2D-image fits.")
        else:
            raise TypeError("Tex should be given as a single value, or 2D numpy.array or 2D-image fits.")

        ##----------------------------
        ##----- loading datacube -----
        if fname[-5:]=='.fits':
            InputMap = fits.open(fname)[0]
            if InputMap.header['NAXIS'] == 4:
                img = InputMap.data[0]
            elif InputMap.header['NAXIS'] == 3:
                img = InputMap.data
            elif (InputMap.header['NAXIS'] < 3) or (InputMap.header['NAXIS'] > 4):
                raise TypeError("data should be given as 3D-spectrum cube fits file.")
            else:
                raise TypeError("data should be given as 3D-spectrum cube fits file.")
        else:
            raise TypeError("data should be given as 3D-spectrum cube fits file.")

        if InputMap.header['CUNIT3'] == 'm s-1':
            onevpix = InputMap.header['CDELT3'] * 0.001  # from m/s to km/s
        elif InputMap.header['CUNIT3'] == 'km s-1':
             onevpix = InputMap.header['CDELT3']         # unit is km/s
        else:
            raise TypeError("The unit of third dimenson has to be km/s or m/s.")
        v0 = InputMap.header['CRVAL3'] * 0.001  # from m/s to km/s
        v0pix = InputMap.header['CRPIX3']
        onepix = 3600. * abs(InputMap.header['CDELT1'])  #pixel size in arcsec
        PixelAreaArcsec = (3600. * abs(InputMap.header['CDELT1'])) * (3600. * abs(InputMap.header['CDELT2']))
        bmaj = InputMap.header['BMAJ'] * 3600
        bmin = InputMap.header['BMIN'] * 3600
        print ( "bmaj=%.4f, bmin=%.4f"%(bmaj, bmin) )
        freq = InputMap.header['RESTFRQ'] / 1e9     # GHz
        if InputMap.header['BUNIT'] == 'beam-1 Jy':
            Jyperbeam2K = 1.222e6/(freq**2 * bmaj * bmin)  # freq in GHz, bmaj and bmin in arcsec
            print ( "Jyperbeam_to_K = %.3f"%(Jyperbeam2K) )
        elif InputMap.header['BUNIT'] == 'K':
            Jyperbeam2K = 1
            print("The unit of datacube is K")
        else:
            raise TypeError("The unit of datacube should be 'Jy/beam' or 'K'.")
        vch_start = np.int( (vl - v0) / onevpix + v0pix )
        vch_stop = np.int( (vu - v0) / onevpix + v0pix )
        print ("----------------------")
        print ('starting velocity', vl,'km/s at channel', vch_start)
        print ('ending velocity  ', vu,'km/s at channel', vch_stop)
        print ("----------------------")

        if vch_start > vch_stop:
            vch_start, vch_stop = vch_stop, vch_start

        ##------load the mask -------
        if mask is not None:
            if (isinstance(mask, np.ndarray)) and (mask.ndim == 2) and (mask.dtype=='int'):
                newmask = mask
            elif (isinstance(mask, np.ndarray)) and (mask.ndim == 2) and (mask.dtype=='bool'):
                newmask = mask.astype(int)      #np.array(mask, dtype=int)  # convert the boolean array to int array
            elif mask[-5:]=='.fits':
                Inputmask = fits.open(mask)[0]
                if Inputmask.data.ndim == 2:
                    newmask = Inputmask.data.astype(int)
                # elif (Inputmask.data.ndim == 2) and (Inputmask.data.dtype=='bool'):
                #     newmask = np.array(Inputmask.data, dtype=int)  #Inputmask.data.astype(int)    # convert the boolean array to int array
                else:
                    raise TypeError("mask should be given as 2D numpy.array or 2D-image fits.")
            else:
                raise TypeError("mask shouldisinstance(pb, (int, float)): be given as 2D numpy.array or 2D-image fits of boolean or int type.")
            ##--- mask the data --------
            img = np.where(newmask == 1, img, np.nan)
            self.mask = newmask
        else:
            warnings.warn("!!! No mask is applied to the datacube")
        self.img = img


        ##----- load the primary beam correction factor ------
        if pb is None:
            pbcor = 1
            print("No primary beam correction is applied to the datacube")
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
            if pbcor.shape != img[vch_start,:,:].shape:
                raise TypeError("'pb' shape is not match to moment0 image")

        emptynarray = np.zeros((img.shape[1],img.shape[2]))
        mass = emptynarray
        momentum = emptynarray
        energy = emptynarray
        flux = emptynarray
        mocol = emptynarray
        for i in range(vch_start, vch_stop+1):
            mom0 = ( img[i,:,:] * np.abs(onevpix) )/pbcor
            fchan =  mom0      #np.sum(mom0[~np.isnan(mom0)])     # Jy/beam
            ##-- compute the total column density of molecule
            N_mol = calcol_line( flux=fchan, name=self.line, J=self.J, Tex=self.Tex, Jyperbeam2K=Jyperbeam2K, tau=self.tau)
            mocol = mocol + N_mol.column
            ##-- compute the gass mass of H2
            Mchan = h2mass(N_mol=N_mol.column, D=self.D, onepix=onepix, Xmol=self.Xmol, decl=self.decl)
            mass = mass + Mchan
            vel = v0 + onevpix * (i-v0pix+1)
            velocity = vel - vlsr
            momentum = momentum + Mchan * np.abs(velocity)
            energy = energy + 0.5 * Mchan * (np.abs(velocity))**2 
            flux = flux + fchan
        self.mass = mass
        self.momentum = momentum
        self.energy = energy
        self.flux = flux
        self.mocol = mocol

