import io
import os
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u

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
class headinfo(object):
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

    head2d: the header corresponing to the 2D image

    wcs2d: the wcs format header corresponing to the 2D image
    '''
    def __init__(self, fname):
        InputMap = fits.open(fname)[0]
        self.head3d = InputMap.header
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
        self.head2d = hduhead
        self.wcs2d = wcs.WCS(hduhead)

##--------------------------------------------------------------------------------------
