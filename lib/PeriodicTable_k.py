#! /usr/bin/python
"""
    periodic table  
"""

Elemass = [1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00,
           20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95,
           39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.85, 58.93,
           58.69, 63.55, 65.47, 69.72, 72.64, 74.92, 78.96, 79.90, 83.80,
           85.47, 87.62, 88.91, 91.22, 92.91, 95.90, 98.00, 101.1, 102.9,
           106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9, 131.3, 
           132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0,
           157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0, 175.0, 178.5, 
           180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.5, 204.4,
           207.2, 209.0, ]

Eledict  = { 'H':1,     'He':2,   'Li':3,    'Be':4,   'B':5,     'C':6,     'N':7,     'O':8,
             'F':9,     'Ne':10,  'Na':11,   'Mg':12,  'Al':13,   'Si':14,   'P':15,    'S':16,
	     'Cl':17,   'Ar':18,  'K':19,    'Ca':20,  'Sc':21,   'Ti':22,   'V':23,    'Cr':24,
	     'Mn':25,   'Fe':26,  'Co':27,   'Ni':28,  'Cu':29,   'Zn':30,   'Ga':31,   'Ge':32,
	     'As':33,   'Se':34,  'Br':35,   'Kr':36,  'Rb':37,   'Sr':38,   'Y':39,    'Zr':40,
	     'Nb':41,   'Mo':42,  'Tc':43,   'Ru':44,  'Rh':45,   'Pd':46,   'Ag':47,   'Cd':48,
	     'In':49,   'Sn':50,  'Sb':51,   'Te':52,  'I':53,    'Xe':54,   'Cs':55,   'Ba':56,
	     'La':57,   'Ce':58,  'Pr':59,   'Nd':60,  'Pm':61,   'Sm':62,   'Eu':63,   'Gd':64, 
	     'Tb':65,   'Dy':66,  'Ho':67,   'Er':68,  'Tm':69,   'Yb':70,   'Lu':71,   'Hf':72, 
	     'Ta':73,   'W':74,   'Re':75,   'Os':76,  'Ir':77,   'Pt':78,   'Au':79,   'Hg':80, 
	     'Tl':81,   'Pb':82,  'Bi':83,   'Po':84,  'At':85,   'Rn':86,   'Fr':87,   'Ra':88, 
	     'Ac':89,   'Th':90,  'Pa':91,   'U':92,   'Np':93,   'Pu':94,   'Am':95,   'Cm':96, 
	     'Bk':97,   'Cf':98,  'Es':99,   'Fm':100, 'Md':101,  'No':102,  'Lr':103,  'Rf':104, 
	     'Db':105,  'Sg':106, 'Bh':107,  'Hs':108, 'Mt':109,  'Ds':110,  'Rg':111,  'Cn':112, 
	     'Uut':113, 'Fl':114, 'Uup':115, 'Lv':116, 'Uus':117, 'UUo':118} 

Eletable = [ 'H'  ,   'He' ,  'Li' ,  'Be' ,  'B'  ,  'C'  ,  'N' ,  'O' , 
             'F'  ,   'Ne' ,  'Na' ,  'Mg' ,  'Al' ,  'Si' ,  'P' ,  'S' ,
	     'Cl' ,   'Ar' ,  'K'  ,  'Ca' ,  'Sc' ,  'Ti' ,  'V' ,  'Cr',
	     'Mn' ,   'Fe' ,  'Co' ,  'Ni' ,  'Cu' ,  'Zn' ,  'Ga',  'Ge',
	     'As' ,   'Se' ,  'Br' ,  'Kr' ,  'Rb' ,  'Sr' ,  'Y' ,  'Zr',
	     'Nb' ,   'Mo' ,  'Tc' ,  'Ru' ,  'Rh' ,  'Pd' ,  'Ag',  'Cd',
	     'In' ,   'Sn' ,  'Sb' ,  'Te' ,  'I'  ,  'Xe' ,  'Cs',  'Ba',
	     'La' ,   'Ce' ,  'Pr' ,  'Nd' ,  'Pm' ,  'Sm' ,  'Eu',  'Gd', 
	     'Tb' ,   'Dy' ,  'Ho' ,  'Er' ,  'Tm' ,  'Yb' ,  'Lu',  'Hf', 
	     'Ta' ,   'W'  ,  'Re' ,  'Os' ,  'Ir' ,  'Pt' ,  'Au',  'Hg', 
	     'Tl' ,   'Pb' ,  'Bi' ,  'Po' ,  'At' ,  'Rn' ,  'Fr',  'Ra', 
	     'Ac' ,   'Th' ,  'Pa' ,  'U'  ,  'Np' ,  'Pu' ,  'Am',  'Cm', 
	     'Bk' ,   'Cf' ,  'Es' ,  'Fm' ,  'Md' ,  'No' ,  'Lr',  'Rf', 
	     'Db' ,   'Sg' ,  'Bh' ,  'Hs' ,  'Mt' ,  'Ds' ,  'Rg',  'Cn', 
	     'Uut',   'Fl' ,  'Uup',  'Lv' ,  'Uus',  'UUo', ] 


if __name__ == '__main__':
    print(Elemass)
