import io
import os
import numpy as np
##--------------------------------------------------------------------------------------
##  Line information
##--------------------------------------------------------------------------------------

line_CO = {'Ju': [1, 2, 3 ,4, 5], 
           'freq': [115.27120, 230.538, 345.79599, 461.04077, 576.26793],  # GHz
           'Eu': [5.5347, 16.59608, 33.1994, 55.3276, 82.9870], #K
           'gu': [3, 5, 7, 9, 11], 
           'Sijmu2':[0.01212, 0.02423, 0.03631, 0.04834, 0.06031],  #D^2
           'B0':[57635.968], # rigrid rotor rotattion constant,# MHz  
           'Ri':[1,1,1,1,1],  
           'note':'diatomic linear molecule'
          }

line_13CO = {'Ju': [1, 2, 3 ], 
           'freq': [110.20135430, 220.39868420, 330.58796530],np
           'Eu': [5.2888, 15.86618, 31.7319], 
           'gu': [3, 5, 7], 
           'Sijmu2':[0.02436, 0.04869, 0.07297], 
           'B0':[55101.01],    
           'Ri': [1,1,1],   
            'note': 'diatomic linear molecule'
          }

line_C18O = {'Ju': [1, 2, 3], 
           'freq': [109.782182, 219.5603541, 329.33055250],
           'Eu': [5.26868, 15.8058, 31.61116], 
           'gu': [3, 5, 7],
           'Sijmu2':[0.01221, 0.0244, 0.03656], 
           'B0':[54891.42],
           'Ri': [1,1,1],
           'note': 'diatomic linear molecule'
          }

line_C17O = {'Ju': [1, 2, 3], 
           'freq': [112.360005, 224.71531, 337.061104],
           'Eu': [5.39234, 16.17688, 32.35322], 
           'gu': [3, 5, 7],
           'Sijmu2':[0.02434, 0.01390, 0.01411],
           'B0':[56179.9911],
           'Ri': [1, 1, 1],
            'note': 'using Ri=1, so need to sum up all hyperfine lines'
          }

line_SiO = {'Ju': [1, 2, 3, 4, 5, 6], 
           'freq': [43.423853, 86.846985, 130.268683, 173.688238, 217.104919, 260.518009],
           'Eu': [2.084, 6.25203, 12.5039, 20.83958, 31.26, 43.76165], 
           'gu': [3, 5, 7, 9, 11, 13], 
           'Sijmu2':[19.26359, 28.89380, 38.52084, 48.14651, 57.77112],
           'B0':[21711.967],
           'Ri': [1,1,1,1,1,1],
           'note': 'diatomic linear molecule'
          }

line_H13COp = {'Ju': [1, 2, 3, 4], 
           'freq': [86.754288, 173.506697, 260.255339, 346.99834400],
           'Eu': [4.16353, 12.49076, 24.98074, 41.63396], 
           'gu': [3, 5, 7, 9], 
           'Sijmu2':[15.21089, 30.41828, 45.62663, 60.84617],
           'B0':[43377.32],
           'Ri': [1,1,1,1,1], ## the hyperfine lines very close, they can be considered as 1 transition
           'note': 'diatomic linear molecule'
          }

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
