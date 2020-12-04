from __future__ import division
import io
import os
import warnings
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u

# import local things
from .calcol import calcol_line

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
    fname:  str
        the name of a 3D-spectrum cube fits file
    
    mask: 2D numpy.array or 2D-image fits of boolean or int type; optional

    pb: a single int or float value, 2D numpy.array, or 2D or 3D fits; optional
        the primary beam correction factor
    
    D: float
        distance in kpc

    Tex: a float value, 2D numpy.array or 2D fits
        the excitation temperature
    
    vlsr: float
        the systematic velocity

    vl, vu: float
        starting and ending velocity
    
    line: str
        the name of molecule
    
    J: str
        the Quantum Numbers (transition) of line; e.g., '1-0', '2-1'
    
    Xmol: float
        the fractional abundance of molecule relative to H2
    
    decl: float, optional
        source declination in degree; default is 0
    
    tau: float, 2D numpy.array or 2D fits; optional
        optical depth; the tau will be used to correct the optical 
        depth effect in colum density if tau is not zero. Defalut is 0 (optical thin).

    -------------- output ----------
    mass: the gas mass

    momentum: the momentum

    energy: the energy

    flux: the integrated flux

    mocol: the molecular column density

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
            raise TypeError("rotational quantum number of molecule should be provided.")

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

        if (InputMap.header['CUNIT3'] == 'm s-1') or (InputMap.header['CUNIT3'] == 'm/s'):
            onevpix = InputMap.header['CDELT3'] * 0.001  # from m/s to km/s
            v0 = InputMap.header['CRVAL3'] * 0.001  # from m/s to km/s
        elif (InputMap.header['CUNIT3'] == 'km s-1') or (InputMap.header['CUNIT3'] == 'km/s'):
             onevpix = InputMap.header['CDELT3']         # unit is km/s
             v0 = InputMap.header['CRVAL3']  # from m/s to km/s
        else:
            raise TypeError("The unit of third dimenson has to be km/s or m/s.")
        
        v0pix = InputMap.header['CRPIX3']
        onepix = 3600. * abs(InputMap.header['CDELT1'])  #pixel size in arcsec
        PixelAreaArcsec = (3600. * abs(InputMap.header['CDELT1'])) * (3600. * abs(InputMap.header['CDELT2']))
        bmaj = InputMap.header['BMAJ'] * 3600
        bmin = InputMap.header['BMIN'] * 3600
        print ( "bmaj=%.4f, bmin=%.4f"%(bmaj, bmin) )
        freq = InputMap.header['RESTFRQ'] / 1e9     # GHz
        if (InputMap.header['BUNIT'] == 'beam-1 Jy') or (InputMap.header['BUNIT'] == 'Jy/beam'):
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
                raise TypeError("mask shoul be given as 2D numpy.array or 2D-image fits of boolean or int type.")
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

