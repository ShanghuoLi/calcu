import io
import os
import numpy as np 

##-----------------------------------------------------------------------------------
#---------------------------
#  Line information
#__________________________
#
# 1 Debye = 1e-18 esu cm  


##------------------------------------------
##-----------------
## 'QNs': rotational quantum number of the lower state in the observed transition
##  Ju=2 for transition J = 2-1, 
## 'freq':line frequency (GHz); 
## 'Eu': line upper energy (Eu/k)
## 'gu': upper state degeneracy 
## 'Aul': Einstein coefficient
## 'B0': rigrid rotor rotattion constant (MHz)
## 'Ri': relative strength R_{i}, The sume of the relative intenisties $\sum_{i}R_{i} = 1$ for each $\Delta J = 1$ transition.
## 'Sijmu2': the product of the transition line strength and the square of the electric dipole ment (debye^2 or just D^2).
## '': 
##-----------------
##------------------------------------------
line_CO = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4'], 
           'freq': [115.27120, 230.538, 345.79599, 461.04077, 576.26793],  # GHz
           'Eu': [5.5347, 16.59608, 33.1994, 55.3276, 82.9870], #K
           'gu': [3, 5, 7, 9, 11], 
          #  'Aul': [10**-7.1425, 10**-6.1605, 10**-5.6026, 10**-5.2128, 10**-4.9132],
          #  'mu': [1.1011e-19],          # molecuar dipole moment, #esu cm
           'Sijmu2':[0.01212, 0.02423, 0.03631, 0.04834, 0.06031],  #D^2
           'B0':[57635.968], # rigrid rotor rotattion constant,# MHz  
           'Ri':[1,1,1,1,1],  
           #relative strength Ri, The sume of the relative intenisties Ri = 1$ for each  J = 1 transition.
           'note':'diatomic linear molecule'
          }
# mu_CO = 1.1011e-19
# S_mu2_CO = (1.1011e-19)**2 = 0.01212 #D^2
# # rigrid rotor rotattion constant,# Hz  
# B0_CO = 57635.968  
# #relative strength R_{i}, The sume of the relative intenisties $\sum_{i}R_{i} = 1$ for each $\Delta J = 1$ transition.
# Ri_CO = 1

# const0 = k/(h*B0_CO)
# # print ( "Qrot = %s*T_ex + 1/3"%(const0) )
# # J = 1
# # K1 = 3 * k / (8 *np.pi**3 * mu_CO**2 * line_CO['QNs'][J] *B0_CO) *10**5 * (k/(h*line_CO['freq'][J]))
# # print ( 'K1=%e cm^-2'%(K1) )

line_13CO = {'QNs': ['1-0', '2-1', '3-2'], 
           'freq': [110.20135430, 220.39868420, 330.58796530],
           'Eu': [5.2888, 15.86618, 31.7319], 
           'gu': [3, 5, 7], 
           'Sijmu2':[0.02436, 0.04869, 0.07297],  #D^2
           'B0':[55101.01],    # rigrid rotor rotattion constant,# MHz  
           'Ri': [1,1,1],    #relative strength R_{i}, The sume of the relative intenisties $\sum_{i}R_{i} = 1$ for each $\Delta J = 1$ transition.
        #    'Aul': [10**-7.1425, 10**-6.1605, 10**-5.6026, 10**-5.2128, 10**-4.9132]  
            'note': 'diatomic linear molecule'
          }


line_C18O = {'QNs': ['1-0', '2-1', '3-2'], 
           'freq': [109.782182, 219.5603541, 329.33055250],
           'Eu': [5.26868, 15.8058, 31.61116], 
           'gu': [3, 5, 7],
           'Sijmu2':[0.01221, 0.0244, 0.03656],  #D^2
           'B0':[54891.42],
           'Ri': [1,1,1],
           'note': 'diatomic linear molecule'
        #    'Aul': [10**-7.1425, 10**-6.1605, 10**-5.6026, 10**-5.2128, 10**-4.9132]  
          }

##-------------------------------
## hyperfine transition
## relative intensities for the correpsonding frequency (transition), 
# if all of the hyperfine were chosen for the calculation Ri=1
##-------------------------------
line_C17O = {'QNs': ['1-0', '2-1', '3-2'], 
           'freq': [112.360005, 224.71531, 337.061104],
           'Eu': [5.39234, 16.17688, 32.35322], 
           'gu': [3, 5, 7],
           'Sijmu2':[0.02434, 0.01390, 0.01411],
           'B0':[56179.9911],
           'Ri': [1, 1, 1],
           ## relative intensities for the correpsonding frequency (transition), 
           # if all of the hyperfine were chosen for the calculation Ri=1
        #    'Aul': [10**-7.1425, 10**-6.1605, 10**-5.6026, 10**-5.2128, 10**-4.9132]  
            'note': 'using Ri=1, so need to sum up all hyperfine lines'
          }


##-------------------------------
## hyperfine transition
## relative intensities for the correpsonding frequency (transition), 
# if all of the hyperfine were chosen for the calculation Ri=1
##-------------------------------
line_13C17O = {'QNs': ['1-0', '2-1', '3-2'], 
           'freq': [107.28855, 214.573873, 321.8515],
           'Eu': [5.14903, 15.44693, 30.89324], 
           'gu': [8, 10, 12],
           'Sijmu2':[0.0326287, 0.04891, 0.06282],
           'B0':[53644.79],
           'Ri': [1, 1, 1],
           ## relative intensities for the correpsonding frequency (transition), 
           # if all of the hyperfine were chosen for the calculation Ri=1
        #    'Aul': [10**-7.1425, 10**-6.1605, 10**-5.6026, 10**-5.2128, 10**-4.9132]  
            'note': 'using Ri=1, so need to sum up all hyperfine lines'
          }


##-------- SiO ---------
line_SiO = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4', '6-5'], 
           'freq': [43.423853, 86.846985, 130.268683, 173.688238, 217.104919, 260.518009],
           'Eu': [2.084, 6.25203, 12.5039, 20.83958, 31.26, 43.76165], 
           'gu': [3, 5, 7, 9, 11, 13], 
           'Sijmu2':[19.26359, 28.89380, 38.52084, 48.14651, 57.77112],
           'B0':[21711.967],
           'Ri': [1,1,1,1,1,1],
           'note': 'diatomic linear molecule'
#            'Aul': [10**-7.1425, 10**-6.1605, 10**-5.6026, 10**-5.2128, 10**-4.9132]  
          }
# # molecuar dipole moment, #esu cm 
# mu_SiO = 3.1e-19  
# S_mu2_SiO = (3.1e-19)**2 
# # rigrid rotor rotattion constant,# Hz  
# B0_SiO = 21711.967  
# #relative strength R_{i}, The sume of the relative intenisties $\sum_{i}R_{i} = 1$ for each $\Delta J = 1$ transition.
# Ri_SiO = 1
# const0 = k/(h*B0_SiO)
# print ( "Qrot = %s*T_ex + 1/3"%(const0) )


## 
line_CS = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4', '6-5', '7-6', '8-7'], 
           'freq': [48.990955, 97.980953,  146.969029, 195.954217, 244.935556, 
            293.912091, 342.882857, 489.750921],
           'Eu': [2.35118, 7.05327, 14.10620, 23.5111, 35.2605, 49.37157, 65.82732, 84.63284], 
          #  'gu': [3, 5, 7, 9, 11, 13, 15, 17], 
           'Sijmu2':[3.82195, 7.05355, 14.10692, 23.51110, 35.26605, 49.37157, 65.82732, 84.63284],
           'B0':[24495.56],
           'Ri': [1,1,1,1,1,1,1,1], 
           'note': ''
          }


line_C34S = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4', '6-5', '7-6', '8-7'], 
           'freq': [48.20694110, 96.41294950, 144.61710070, 192.81845660, 241.01608920, 289.20906840, 
           337.39645900, 385.57733600],
           'Eu': [2.31356, 6.24804, 11.80342, 18.97955, 27.77627, 38.19355, 50.23093, 63.88848], 
          #  'gu': [6, 10, 14, 18, 22, 26, 30, 34], 
           'Sijmu2':[3.83344, 7.64989, 11.42218, 15.12507, 18.73067, 22.21857, 25.56800, 28.74990],
           'B0':[24103.54],
           'Ri': [1,1,1,1,1,1,1,1], 
           'note': ''
          }

## hyperfine transition
line_N2Hp = {'QNs': ['1-0'], 
           'freq': [93.1733977],
           'Eu': [4.4716], 
           'gu': [3], 
           'Sijmu2':[104.03991770513517],
           'B0':[46586.87],
           'Ri': [5.0/9], 
           'note': 'using the mian hyperfine line, Ri=5/9; only for J=1-0'
          }

## hyperfine transition
line_N2Dp = {'QNs': ['1-0', '2-1', '3-2' ,'4-3'], 
           'freq': [77.1092433, 154.2170112, 231.3218283, 385.5167210],
           'Eu': [3.701, 11.10177, 22.20337, 37.00526], 
           'gu': [3, 5, 7, 9], 
           'Sijmu2':[104.03992, 208.06546, 312.10448, 520.2355186873041],
           'B0':[38554.74],
           'Ri': [0,0,1,1], 
           'note': '''For J=1-0 and 2-1, 
               the separation between the hyperfine lines is several km/s, 
               therefore the Ri is depending on which hyperfine line is used. 
               It is Ri=1 if sum up all hyperfine lines. 
                You can take a look at Redaelli, E., et al 2019'''
          }

## hyperfine transition
line_DCN = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4'], 
           'freq': [72.4146936, 144.8280015, 217.2385378, 289.6449170, 362.0457535],
           'Eu': [3.47528, 10.42595, 20.85164, 34.75234, 52.12784], 
           'gu': [3, 5, 7, 9, 11], 
           'Sijmu2':[26.83490115967991, 53.678454102551235, 80.50076361080625, 107.3494935805236, 134.1837732253372],
           'B0':[36207.46],
           'Ri': [1,1,1,1,1], 
           'note': ''
          }

## hyperfine transition
line_DCOp = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4'], 
           'freq': [72.0393124, 144.0772890, 216.1125822, 288.1438583, 360.1697783],
           'Eu': [3.457, 10.37, 20.74, 34.57224, 51.85764], 
           'gu': [3, 5, 7, 9, 11], 
           'Sijmu2':[15.212479290080404, 30.42285727480961, 45.62469770768503, 60.84151222917087, 76.05013686927006],
           'B0':[46586.88],
           'Ri': [1,1,1,1,1], 
           'note': ''
          }

## hyperfine transition
line_13CS = {'QNs': ['1-0', '2-1', '3-2' ,'4-3', '5-4', '6-5', '7-6', '8-7'], 
           'freq': [46.2475632, 92.4943080, 138.7393350, 184.9817720, 231.2206852, 
           277.4554050, 323.6849730, 369.9085505],
           'Eu': [2.21952, 6.65859,13.31687,22.19462,33.29137,46.60706,62.14138,79.89413], 
           'gu': [6, 10, 14, 18, 22, 26, 30, 34], 
           'Sijmu2':[7.66726, 15.33557, 23.00382, 30.669023, 
           38.33540, 46.00447, 53.666096, 61.34794],
           'B0':[46586.88],
           'Ri': [1,1,1,1,1,1,1,1], 
           'note': ''
          }

## hyperfine transition
line_C2D = {'QNs': ['1-0', '2-1', '3-2'], 
           'freq': [0, 0, 216.37283],
           'Eu': [3.345, 10.383, 20.77], 
           'gu': [3, 5, 7, 9, 11], 
           'Sijmu2':[0, 0, 2.54121],
           'B0':[36068.02],
           'Ri': [1,1,1], 
           'note': 'only J=3-2 available'
          }

#
#
##-------- H13CO+ ---------
line_H13COp = {'QNs': ['1-0', '2-1', '3-2' ,'4-3'], 
           'freq': [86.754288, 173.506697, 260.255339, 346.99834400],
           'Eu': [4.16353, 12.49076, 24.98074, 41.63396], 
           'gu': [3, 5, 7, 9], 
           'Sijmu2':[15.21089, 30.41828, 45.62663, 60.84617],
           'B0':[43377.32],
           'Ri': [1,1,1,1,1], ## the hyperfine lines very close, they can be considered as 1 transition
           'note': 'diatomic linear molecule'
#            'Aul': []  
          }
# molecuar dipole moment, #esu cm 
# mu_H13COP = 3.888e-19
# S_mu2_H13COP = 15.21089 #D^2
# # rigrid rotor rotattion constant,# Hz  
# B0_H13COP = 43377.32  
# #relative strength R_{i}, The sume of the relative intenisties $\sum_{i}R_{i} = 1$ for each $\Delta J = 1$ transition.
# Ri_H13COP = 1self.molecule['B0'][0]



##-------- NH2D ---------
line_NH2D = {'QNs': ['1(1,1)-1(0,1)'], 
           'freq': [85.926],
           'Eu': [20.68], 
           'gu': [15], 
           'Sijmu2': [11.91],
           'Ri': [1/2], #the relative intensity of the main hyperfine transition with respect to the others.  
           'note': '''using the main hyperfine only, the relative intensity of the main 
           hyperfine transition with respect to the others'''
          }
# molecuar dipole moment, #esu cm 
# mu_NH2D = 3.888e-19
# S_mu2_NH2D = 11.91
# # rigrid rotor rotattion constant,# Hz  
# # B0_NH2D = 140795
# #relative strength R_{i}, The sume of the relative intenisties $\sum_{i}R_{i} = 1$ for each $\Delta J = 1$ transition.
# Ri_NH2D = 1/2



##-----------------------------------------------------------------------
## slightly-asymmetric rotor molecule
##-----------------------------------------------------------------------
# 
line_H2CO = {'QNs': ['3(0,3)-2(0,2)', '3(2,2)-2(2,1)', '3(2,1)-2(2,0)',
                    '4(0,4)-3(0,3)', '4(2,3)-3(2,2)', '4(2,2)-3(2,1)'], 
           'freq': [218.222192, 218.4756320, 218.7600660,
                    290.62341, 291.23778, 291.94806],
           'Eu': [20.9564, 68.0937, 68.11081, 34.90, 82.07, 82.12], 
           'gu': [7, 7, 7], 
           'Sijmu2':[16.307973610235937, 9.062123990047878,9.062123990047878,
                      21.74202546161732, 16.307973610235937,16.31172909227838],
           'B0':[281970.56,38833.987,34004.244],
           'Ri': [1,1,1, 1,1,1], 
          #  'Ri': [3./7, 5./21, 5./21, 4./9, 1./3, 1./3], 
           'note': ''
          }
##3(0,3)-2(0,2): J(K_A, K_C)-(J-1)(K_A, K_C-1)

#
line_CH3OH = {'QNs': ['4(-2,3)-3(-1,2)E'], 
           'freq': [218.440063],
           'Eu': [45.45944], 
           'Sijmu2':[13.905928754654452],
           'B0':[127523.4, 24690.2, 23759.7],
           'Ri': [1], 
           'note': ''
          }


line_HNCO = {'QNs': ['10(0,10)-9(0,9)'], 
           'freq': [219.798274],
           'Eu': [58.], 
           'Sijmu2':[24.963185272433364 ],
           'B0':[912711.4, 11071.0, 10910.57],
           'Ri': [1], 
           'note': ''
          }

