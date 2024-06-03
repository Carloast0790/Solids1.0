# get_chemical_symbol(sym/z) : RETURNS A SYMBOL, BASED ON THE DICTIONARY chemical_symbol[z]
# get_atomic_mass(sym)       : RETURNS THE ATOMIC MASS IN AMU, BASED ON THE atomic_mass[sym]
# get_covalent_radius(sym)   : RETURNS THE COVALENT RADIUS IN A, BASED ON THE covalent_radius[sym]
# check_atomic_sym(sym)      : RETURNS TRUE OF FALSE, BASED ON THE chemical_symbols
#------------------------------------------------------------------------------------------
def is_number(s):
    try:
        float(s).is_integer()
        return True
    except ValueError:
        pass
#------------------------------------------------------------------------------------------
""" 
MAXIMUM ATOMIC NUMBER (Z) IS 86
ATOMIC SYMBOL IS USED BY DEFAULT
ROBUSTS, sym COULD BE a SYMBOL ("H") OR an INTEGER (1)
HOW TO USE:
>>> from utils.atomic import get_chemical_symbol
"""

def get_chemical_symbol(sym):
    chemical_symbol={}
    chemical_symbol[1]="H"
    chemical_symbol[2]="He"
    chemical_symbol[3]="Li"
    chemical_symbol[4]="Be"
    chemical_symbol[5]="B"
    chemical_symbol[6]="C"
    chemical_symbol[7]="N"
    chemical_symbol[8]="O"
    chemical_symbol[9]="F"
    chemical_symbol[10]="Ne"
    chemical_symbol[11]="Na"
    chemical_symbol[12]="Mg"
    chemical_symbol[13]="Al"
    chemical_symbol[14]="Si"
    chemical_symbol[15]="P"
    chemical_symbol[16]="S"
    chemical_symbol[17]="Cl"
    chemical_symbol[18]="Ar"
    chemical_symbol[19]="K"
    chemical_symbol[20]="Ca"
    chemical_symbol[21]="Sc"
    chemical_symbol[22]="Ti"
    chemical_symbol[23]="V"
    chemical_symbol[24]="Cr"
    chemical_symbol[25]="Mn"
    chemical_symbol[26]="Fe"
    chemical_symbol[27]="Co"
    chemical_symbol[28]="Ni"
    chemical_symbol[29]="Cu"
    chemical_symbol[30]="Zn"
    chemical_symbol[31]="Ga"
    chemical_symbol[32]="Ge"
    chemical_symbol[33]="As"
    chemical_symbol[34]="Se"
    chemical_symbol[35]="Br"
    chemical_symbol[36]="Kr"
    chemical_symbol[37]="Rb"
    chemical_symbol[38]="Sr"
    chemical_symbol[39]="Y"
    chemical_symbol[40]="Zr"
    chemical_symbol[41]="Nb"
    chemical_symbol[42]="Mo"
    chemical_symbol[43]="Tc"
    chemical_symbol[44]="Ru"
    chemical_symbol[45]="Rh"
    chemical_symbol[46]="Pd"
    chemical_symbol[47]="Ag"
    chemical_symbol[48]="Cd"
    chemical_symbol[49]="In"
    chemical_symbol[50]="Sn"
    chemical_symbol[51]="Sb"
    chemical_symbol[52]="Te"
    chemical_symbol[53]="I"
    chemical_symbol[54]="Xe"
    chemical_symbol[55]="Cs"
    chemical_symbol[56]="Ba"
    chemical_symbol[57]="La"
    chemical_symbol[58]="Ce"
    chemical_symbol[59]="Pr"
    chemical_symbol[60]="Nd"
    chemical_symbol[61]="Pm"
    chemical_symbol[62]="Sm"
    chemical_symbol[63]="Eu"
    chemical_symbol[64]="Gd"
    chemical_symbol[65]="Tb"
    chemical_symbol[66]="Dy"
    chemical_symbol[67]="Ho"
    chemical_symbol[68]="Er"
    chemical_symbol[69]="Tm"
    chemical_symbol[70]="Yb"
    chemical_symbol[71]="Lu"
    chemical_symbol[72]="Hf"
    chemical_symbol[73]="Ta"
    chemical_symbol[74]="W"
    chemical_symbol[75]="Re"
    chemical_symbol[76]="Os"
    chemical_symbol[77]="Ir"
    chemical_symbol[78]="Pt"
    chemical_symbol[79]="Au"
    chemical_symbol[80]="Hg"
    chemical_symbol[81]="Tl"
    chemical_symbol[82]="Pb"
    chemical_symbol[83]="Bi"
    chemical_symbol[84]="Po"
    chemical_symbol[85]="At"
    chemical_symbol[86]="Rn"
    chemical_symbol[87]="Fr"
    chemical_symbol[88]="Ra"
    chemical_symbol[89]="Ac"
    chemical_symbol[90]="Th"
    chemical_symbol[91]="Pa"
    chemical_symbol[92]="U"
    chemical_symbol[93]="Np"
    chemical_symbol[94]="Pu"
    chemical_symbol[95]="Am"
    chemical_symbol[96]="Cm"
    chemical_symbol[97]="Bk"
    chemical_symbol[98]="Cf"
    chemical_symbol[99]="Es"
    chemical_symbol[100]="Fm"
    chemical_symbol[101]="Md"
    chemical_symbol[102]="No"
    chemical_symbol[103]="Lr"
    chemical_symbol[104]="Rf"
    chemical_symbol[105]="Db"
    chemical_symbol[106]="Sg"
    chemical_symbol[107]="Bh"
    chemical_symbol[108]="Hs"
    chemical_symbol[109]="Mt"
    chemical_symbol[110]="Ds"
    chemical_symbol[111]="Rg"
    chemical_symbol[112]="Cn"
    chemical_symbol[113]="Nh"
    chemical_symbol[114]="Fl"
    chemical_symbol[115]="Mc"
    chemical_symbol[116]="Lv"
    chemical_symbol[117]="Ts"
    chemical_symbol[118]="Og"
    sym1 = chemical_symbol[int(sym)] if is_number(sym) else sym
    return sym1
#------------------------------------------------------------------------------------------
"""
MAXIMUM ATOMIC NUMBER (Z) IS 86
ATOMIC SYMBOL IS USED BY DEFAULT
HOW TO USE:
>>> from utils.atomic import check_atomic_sym
"""

def check_atomic_sym(sym):
    val0=False
    for ii in range(1, 87):
        if sym == get_chemical_symbol(ii): val0=True
    return val0
#------------------------------------------------------------------------------------------    
"""
MAXIMUM ATOMIC NUMBER (Z) IS 86
ATOMIC SYMBOL IS USED BY DEFAULT
HOW TO USE:
>>> from utils.atomic import get_atomic_number
"""

def get_atomic_number(sym):
    for ii in range(1, 87):
        if sym == get_chemical_symbol(ii): atomicnumber=int(ii)
    return atomicnumber
#------------------------------------------------------------------------------------------
"""
MAXIMUM ATOMIC NUMBER (Z) IS 86
ATOMIC SYMBOL IS USED BY DEFAULT
Units: g/mol = amu
HOW TO USE:
>>> from utils.atomic import get_atomic_mass
https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
"""

def get_atomic_mass(sym):
    relative_atomic_mass={}
    relative_atomic_mass["H"]=1.00782503223
    relative_atomic_mass["He"]=3.0160293201
    relative_atomic_mass["Li"]=6.0151228874
    relative_atomic_mass["Be"]=9.012183065
    relative_atomic_mass["B"]=10.01293695
    relative_atomic_mass["C"]=12.0000000
    relative_atomic_mass["N"]=14.00307400443
    relative_atomic_mass["O"]=15.99491461957
    relative_atomic_mass["F"]=18.99840316273
    relative_atomic_mass["Ne"]=19.9924401762
    relative_atomic_mass["Na"]=22.9897692820
    relative_atomic_mass["Mg"]=23.985041697
    relative_atomic_mass["Al"]=26.98153853
    relative_atomic_mass["Si"]=27.97692653465
    relative_atomic_mass["P"]=30.97376199842
    relative_atomic_mass["S"]=31.9720711744
    relative_atomic_mass["Cl"]=34.968852682
    relative_atomic_mass["Ar"]=35.967545105
    relative_atomic_mass["K"]=38.9637064864
    relative_atomic_mass["Ca"]=39.962590863
    relative_atomic_mass["Sc"]=44.95590828
    relative_atomic_mass["Ti"]=45.95262772
    relative_atomic_mass["V"]=49.94715601
    relative_atomic_mass["Cr"]=49.94604183
    relative_atomic_mass["Mn"]=54.93804391
    relative_atomic_mass["Fe"]=53.93960899
    relative_atomic_mass["Co"]=58.93319429
    relative_atomic_mass["Ni"]=57.93534241
    relative_atomic_mass["Cu"]=62.92959772
    relative_atomic_mass["Zn"]=63.92914201
    relative_atomic_mass["Ga"]=68.9255735
    relative_atomic_mass["Ge"]=69.92424875
    relative_atomic_mass["As"]=74.92159457
    relative_atomic_mass["Se"]=73.922475934
    relative_atomic_mass["Br"]=78.9183376
    relative_atomic_mass["Kr"]=77.92036494
    relative_atomic_mass["Rb"]=84.9117897379
    relative_atomic_mass["Sr"]=83.9134191
    relative_atomic_mass["Y"]=88.9058403
    relative_atomic_mass["Zr"]=89.9046977
    relative_atomic_mass["Nb"]=92.9063730
    relative_atomic_mass["Mo"]=91.90680796
    relative_atomic_mass["Tc"]=96.9063667
    relative_atomic_mass["Ru"]=95.90759025
    relative_atomic_mass["Rh"]=102.9054980
    relative_atomic_mass["Pd"]=101.9056022
    relative_atomic_mass["Ag"]=106.9050916
    relative_atomic_mass["Cd"]=105.9064599
    relative_atomic_mass["In"]=112.90406184
    relative_atomic_mass["Sn"]=111.90482387
    relative_atomic_mass["Sb"]=120.9038120
    relative_atomic_mass["Te"]=119.9040593
    relative_atomic_mass["I"]=126.9044719
    relative_atomic_mass["Xe"]=123.9058920
    relative_atomic_mass["Cs"]=132.9054519610
    relative_atomic_mass["Ba"]=129.9063207
    relative_atomic_mass["La"]=137.9071149
    relative_atomic_mass["Ce"]=135.90712921
    relative_atomic_mass["Pr"]=140.9076576
    relative_atomic_mass["Nd"]=141.9077290
    relative_atomic_mass["Pm"]=144.9127559
    relative_atomic_mass["Sm"]=143.9120065
    relative_atomic_mass["Eu"]=150.9198578
    relative_atomic_mass["Gd"]=151.9197995
    relative_atomic_mass["Tb"]=158.9253547
    relative_atomic_mass["Dy"]=155.9242847
    relative_atomic_mass["Ho"]=164.9303288
    relative_atomic_mass["Er"]=161.9287884
    relative_atomic_mass["Tm"]=168.9342179
    relative_atomic_mass["Yb"]=167.9338896
    relative_atomic_mass["Lu"]=174.9407752
    relative_atomic_mass["Hf"]=173.9400461
    relative_atomic_mass["Ta"]=179.9474648
    relative_atomic_mass["W"]=179.9467108
    relative_atomic_mass["Re"]=184.9529545
    relative_atomic_mass["Os"]=183.9524885
    relative_atomic_mass["Ir"]=190.9605893
    relative_atomic_mass["Pt"]=189.9599297
    relative_atomic_mass["Au"]=196.96656879
    relative_atomic_mass["Hg"]=195.9658326
    relative_atomic_mass["Tl"]=202.9723446
    relative_atomic_mass["Pb"]=203.9730440
    relative_atomic_mass["Bi"]=208.9803991
    relative_atomic_mass["Po"]=208.9824308
    relative_atomic_mass["At"]=209.9871479
    relative_atomic_mass["Rn"]=210.9906011
    relative_atomic_mass["Fr"]=223.0197360
    relative_atomic_mass["Ra"]=223.0185023
    relative_atomic_mass["Ac"]=227.0277523
    relative_atomic_mass["Th"]=230.0331341
    relative_atomic_mass["Pa"]=231.0358842
    relative_atomic_mass["U"]=233.0396355
    relative_atomic_mass["Np"]=236.046570
    relative_atomic_mass["Pu"]=238.0495601
    relative_atomic_mass["Am"]=241.0568293
    relative_atomic_mass["Cm"]=243.0613893
    relative_atomic_mass["Bk"]=247.0703073
    relative_atomic_mass["Cf"]=249.0748539
    relative_atomic_mass["Es"]=252.082980
    relative_atomic_mass["Fm"]=257.0951061
    relative_atomic_mass["Md"]=258.0984315
    relative_atomic_mass["No"]=259.10103
    relative_atomic_mass["Lr"]=262.10961
    relative_atomic_mass["Rf"]=267.12179
    relative_atomic_mass["Db"]=268.12567
    relative_atomic_mass["Sg"]=271.13393
    relative_atomic_mass["Bh"]=272.13826
    relative_atomic_mass["Hs"]=270.13429
    relative_atomic_mass["Mt"]=276.15159
    relative_atomic_mass["Ds"]=281.16451
    relative_atomic_mass["Rg"]=280.16514
    relative_atomic_mass["Cn"]=285.17712
    relative_atomic_mass["Nh"]=284.17873
    relative_atomic_mass["Fl"]=289.19042
    relative_atomic_mass["Mc"]=288.19274
    relative_atomic_mass["Lv"]=293.20449
    relative_atomic_mass["Ts"]=292.20746
    relative_atomic_mass["Og"]=294.21392
    atomicmass=relative_atomic_mass[sym]
    return atomicmass
#------------------------------------------------------------------------------------------
def get_valence_electrons(sym):
    valence_electrons={}
    valence_electrons["H"] = 1.0
    valence_electrons["He"]= 8.0
    valence_electrons["Li"]= 1.0
    valence_electrons["Be"]= 2.0
    valence_electrons["B"] = 3.0
    valence_electrons["C"] = 4.0
    valence_electrons["N"] = 5.0
    valence_electrons["O"] = 6.0
    valence_electrons["F"] = 7.0
    valence_electrons["Ne"]= 10.0
    valence_electrons["Na"]= 1.0
    valence_electrons["Mg"]= 2.0
    valence_electrons["Al"]= 3.0
    valence_electrons["Si"]= 4.0
    valence_electrons["P"] = 5.0
    valence_electrons["S"] = 6.0
    valence_electrons["Cl"]= 7.0
    valence_electrons["Ar"]= 8.0
    valence_electrons["K"] = 1.0
    valence_electrons["Ca"]= 2.0
    valence_electrons["Sc"]= 3.0
    valence_electrons["Ti"]= 4.0
    valence_electrons["V"] = 5.0
    valence_electrons["Cr"]= 6.0
    valence_electrons["Mn"]= 7.0
    valence_electrons["Fe"]= 8.0
    valence_electrons["Co"]= 8.0
    valence_electrons["Ni"]= 8.0
    valence_electrons["Cu"]= 1.0
    valence_electrons["Zn"]= 2.0
    valence_electrons["Ga"]= 3.0
    valence_electrons["Ge"]= 4.0
    valence_electrons["As"]= 5.0
    valence_electrons["Se"]= 6.0
    valence_electrons["Br"]= 7.0
    valence_electrons["Kr"]= 8.0
    valence_electrons["Rb"]= 1.0
    valence_electrons["Sr"]= 2.0
    valence_electrons["Y"] = 3.0
    valence_electrons["Zr"]= 4.0
    valence_electrons["Nb"]= 5.0
    valence_electrons["Mo"]= 6.0
    valence_electrons["Tc"]= 7.0
    valence_electrons["Ru"]= 8.0
    valence_electrons["Rh"]= 8.0
    valence_electrons["Pd"]= 8.0
    valence_electrons["Ag"]= 1.0
    valence_electrons["Cd"]= 2.0
    valence_electrons["In"]= 3.0
    valence_electrons["Sn"]= 4.0
    valence_electrons["Sb"]= 5.0
    valence_electrons["Te"]= 6.0
    valence_electrons["I"] = 7.0
    valence_electrons["Xe"]= 8.0
    valence_electrons["Cs"]= 1.0
    valence_electrons["Ba"]= 2.0
    valence_electrons["La"]= 3.0
    valence_electrons["Ce"]= 3.0
    valence_electrons["Pr"]= 3.0
    valence_electrons["Nd"]= 3.0
    valence_electrons["Pm"]= 3.0
    valence_electrons["Sm"]= 3.0
    valence_electrons["Eu"]= 3.0
    valence_electrons["Gd"]= 3.0
    valence_electrons["Tb"]= 3.0
    valence_electrons["Dy"]= 3.0
    valence_electrons["Ho"]= 3.0
    valence_electrons["Er"]= 3.0
    valence_electrons["Tm"]= 3.0
    valence_electrons["Yb"]= 3.0
    valence_electrons["Lu"]= 3.0
    valence_electrons["Hf"]= 4.0
    valence_electrons["Ta"]= 5.0
    valence_electrons["W"] = 6.0
    valence_electrons["Re"]= 7.0
    valence_electrons["Os"]= 8.0
    valence_electrons["Ir"]= 8.0
    valence_electrons["Pt"]= 8.0
    valence_electrons["Au"]= 1.0
    valence_electrons["Hg"]= 2.0
    valence_electrons["Tl"]= 3.0
    valence_electrons["Pb"]= 4.0
    valence_electrons["Bi"]= 5.0
    valence_electrons["Po"]= 6.0
    valence_electrons["At"]= 7.0
    valence_electrons["Rn"]= 8.0
    valence_e=valence_electrons[sym]
    return valence_e
#------------------------------------------------------------------------------------------

"""
MAXIMUM ATOMIC NUMBER (Z) IS 86
ATOMIC SYMBOL IS USED BY DEFAULT
Units : Angstrom
HOW TO USE:
>>> from utils.atomic import get_covalent_radius
"""

def get_covalent_radius(sym):
    covalent_radius={}
    covalent_radius["H"] =0.32
    covalent_radius["He"]=0.93
    covalent_radius["Li"]=1.23
    covalent_radius["Be"]=0.90
    covalent_radius["B"] =0.82
    covalent_radius["C"] =0.77
    covalent_radius["N"] =0.75
    covalent_radius["O"] =0.73
    covalent_radius["F"] =0.72
    covalent_radius["Ne"]=0.71
    covalent_radius["Na"]=1.54
    covalent_radius["Mg"]=1.36
    covalent_radius["Al"]=1.18
    covalent_radius["Si"]=1.11
    covalent_radius["P"] =1.06
    covalent_radius["S"] =1.02
    covalent_radius["Cl"]=0.99
    covalent_radius["Ar"]=0.98
    covalent_radius["K"] =2.03
    covalent_radius["Ca"]=1.74
    covalent_radius["Sc"]=1.44
    covalent_radius["Ti"]=1.32
    covalent_radius["V"] =1.22
    covalent_radius["Cr"]=1.18
    covalent_radius["Mn"]=1.17
    covalent_radius["Fe"]=1.17
    covalent_radius["Co"]=1.16
    covalent_radius["Ni"]=1.15
    covalent_radius["Cu"]=1.17
    covalent_radius["Zn"]=1.25
    covalent_radius["Ga"]=1.26
    covalent_radius["Ge"]=1.22
    covalent_radius["As"]=1.20
    covalent_radius["Se"]=1.16
    covalent_radius["Br"]=1.14
    covalent_radius["Kr"]=1.12
    covalent_radius["Rb"]=2.16
    covalent_radius["Sr"]=1.91
    covalent_radius["Y"] =1.62
    covalent_radius["Zr"]=1.45
    covalent_radius["Nb"]=1.34
    covalent_radius["Mo"]=1.30
    covalent_radius["Tc"]=1.27
    covalent_radius["Ru"]=1.25
    covalent_radius["Rh"]=1.25
    covalent_radius["Pd"]=1.28
    covalent_radius["Ag"]=1.34
    covalent_radius["Cd"]=1.48
    covalent_radius["In"]=1.44
    covalent_radius["Sn"]=1.41
    covalent_radius["Sb"]=1.40
    covalent_radius["Te"]=1.36
    covalent_radius["I"] =1.33
    covalent_radius["Xe"]=1.31
    covalent_radius["Cs"]=2.35
    covalent_radius["Ba"]=1.98
    covalent_radius["La"]=1.69
    covalent_radius["Ce"]=1.65
    covalent_radius["Pr"]=1.65
    covalent_radius["Nd"]=1.64
    covalent_radius["Pm"]=1.63
    covalent_radius["Sm"]=1.62
    covalent_radius["Eu"]=1.85
    covalent_radius["Gd"]=1.61
    covalent_radius["Tb"]=1.59
    covalent_radius["Dy"]=1.59
    covalent_radius["Ho"]=1.58
    covalent_radius["Er"]=1.57
    covalent_radius["Tm"]=1.56
    covalent_radius["Yb"]=1.74
    covalent_radius["Lu"]=1.56
    covalent_radius["Hf"]=1.44
    covalent_radius["Ta"]=1.34
    covalent_radius["W"] =1.30
    covalent_radius["Re"]=1.28
    covalent_radius["Os"]=1.26
    covalent_radius["Ir"]=1.27
    covalent_radius["Pt"]=1.30
    covalent_radius["Au"]=1.34
    covalent_radius["Hg"]=1.49
    covalent_radius["Tl"]=1.48
    covalent_radius["Pb"]=1.47
    covalent_radius["Bi"]=1.46
    covalent_radius["Po"]=1.46
    covalent_radius["At"]=1.45
    covalent_radius["Rn"]=1.90
    covalentradius=covalent_radius[sym]
    return covalentradius
#------------------------------------------------------------------------------------------
"""
MAXIMUM ATOMIC NUMBER (Z) IS 86
ATOMIC SYMBOL IS USED BY DEFAULT
HOW TO USE:
>>> from utils.atomic import get_chemical_color
"""

def get_chemical_color(sym):
    chemical_color={}
    chemical_color['H'] ='<1.00000000,1.00000000,1.00000000>'
    chemical_color['He']='<0.85098039,1.00000000,1.00000000>'
    chemical_color['Li']='<0.80000000,0.50196078,1.00000000>'
    chemical_color['Be']='<0.76078431,1.00000000,0.00000000>'
    chemical_color['B'] ='<1.00000000,0.70980392,0.70980392>'
    chemical_color['C'] ='<0.56470588,0.56470588,0.56470588>'
    chemical_color['N'] ='<0.18823529,0.31372549,0.97254902>'
    chemical_color['O'] ='<1.00000000,0.05098039,0.05098039>'
    chemical_color['F'] ='<0.56470588,0.87843137,0.31372549>'
    chemical_color['Ne']='<0.70196078,0.89019608,0.96078431>'
    chemical_color['Na']='<0.67058824,0.36078431,0.94901961>'
    chemical_color['Mg']='<0.54117647,1.00000000,0.00000000>'
    chemical_color['Al']='<0.74901961,0.65098039,0.65098039>'
    chemical_color['Si']='<0.94117647,0.78431373,0.62745098>'
    chemical_color['P'] ='<1.00000000,0.50196078,0.00000000>'
    chemical_color['S'] ='<1.00000000,1.00000000,0.18823529>'
    chemical_color['Cl']='<0.12156863,0.94117647,0.12156863>'
    chemical_color['Ar']='<0.50196078,0.81960784,0.89019608>'
    chemical_color['K'] ='<0.56078431,0.25098039,0.83137255>'
    chemical_color['Ca']='<0.23921569,1.00000000,0.00000000>'
    chemical_color['Sc']='<0.90196078,0.90196078,0.90196078>'
##  chemical_color['Ti']='<0.74901961,0.76078431,0.78039216>'
    chemical_color['Ti']='<0.45100000,0.76100000,0.98400000>'
    chemical_color['V'] ='<0.65098039,0.65098039,0.67058824>'
    chemical_color['Cr']='<0.54117647,0.60000000,0.78039216>'
    chemical_color['Mn']='<0.61176471,0.47843137,0.78039216>'
    chemical_color['Fe']='<0.87843137,0.40000000,0.20000000>'
    chemical_color['Co']='<0.94117647,0.56470588,0.62745098>'
    chemical_color['Ni']='<0.31372549,0.81568627,0.31372549>'
    chemical_color['Cu']='<0.78431373,0.50196078,0.20000000>'
    chemical_color['Zn']='<0.49019608,0.50196078,0.69019608>'
    chemical_color['Ga']='<0.76078431,0.56078431,0.56078431>'
    chemical_color['Ge']='<0.40000000,0.56078431,0.56078431>'
    chemical_color['As']='<0.74117647,0.50196078,0.89019608>'
    chemical_color['Se']='<1.00000000,0.63137255,0.00000000>'
    chemical_color['Br']='<0.65098039,0.16078431,0.16078431>'
    chemical_color['Kr']='<0.36078431,0.72156863,0.81960784>'
    chemical_color['Rb']='<0.43921569,0.18039216,0.69019608>'
    chemical_color['Sr']='<0.00000000,1.00000000,0.00000000>'
    chemical_color['Y'] ='<0.58039216,1.00000000,1.00000000>'
##  chemical_color['Zr']='<0.58039215,0.87843137,0.87843137>'
    chemical_color['Zr']='<0.0,1.0,0.5>'
    chemical_color['Nb']='<0.45098039,0.76078431,0.78823529>'
    chemical_color['Mo']='<0.32941176,0.70980392,0.70980392>'
    chemical_color['Tc']='<0.23137255,0.61960784,0.61960784>'
    chemical_color['Ru']='<0.14117647,0.56078431,0.56078431>'
    chemical_color['Rh']='<0.03921569,0.49019608,0.54901961>'
    chemical_color['Pd']='<0.00000000,0.41176471,0.52156863>'
    chemical_color['Ag']='<0.75294118,0.75294118,0.75294118>'
    chemical_color['Cd']='<1.00000000,0.85098039,0.56078431>'
    chemical_color['In']='<0.65098039,0.45882353,0.45098039>'
    chemical_color['Sn']='<0.40000000,0.50196078,0.50196078>'
    chemical_color['Sb']='<0.61960784,0.38823529,0.70980392>'
    chemical_color['Te']='<0.83137255,0.47843137,0.00000000>'
    chemical_color['I'] ='<0.58039216,0.00000000,0.58039216>'
    chemical_color['Xe']='<0.25882353,0.61960784,0.69019608>'
    chemical_color['Cs']='<0.34117647,0.09019608,0.56078431>'
    chemical_color['Ba']='<0.00000000,0.78823529,0.00000000>'
    chemical_color['La']='<0.43921569,0.83137255,1.00000000>'
    chemical_color['Ce']='<1.00000000,1.00000000,0.78039216>'
    chemical_color['Pr']='<0.85098039,1.00000000,0.78039216>'
    chemical_color['Nd']='<0.78039216,1.00000000,0.78039216>'
    chemical_color['Pm']='<0.63921569,1.00000000,0.78039216>'
    chemical_color['Sm']='<0.56078431,1.00000000,0.78039216>'
    chemical_color['Eu']='<0.38039216,1.00000000,0.78039216>'
    chemical_color['Gd']='<0.27058824,1.00000000,0.78039216>'
    chemical_color['Tb']='<0.18823529,1.00000000,0.78039216>'
    chemical_color['Dy']='<0.12156863,1.00000000,0.78039216>'
    chemical_color['Ho']='<0.00000000,1.00000000,0.61176471>'
    chemical_color['Er']='<0.00000000,0.90196078,0.45882353>'
    chemical_color['Tm']='<0.00000000,0.83137255,0.32156863>'
    chemical_color['Yb']='<0.00000000,0.74901961,0.21960784>'
    chemical_color['Lu']='<0.00000000,0.67058824,0.14117647>'
    chemical_color['Hf']='<0.30196078,0.76078431,1.00000000>'
    chemical_color['Ta']='<0.30196078,0.65098039,1.00000000>'
    chemical_color['W'] ='<0.12941176,0.58039216,0.83921569>'
    chemical_color['Re']='<0.14901961,0.49019608,0.67058824>'
    chemical_color['Os']='<0.14901961,0.40000000,0.58823529>'
    chemical_color['Ir']='<0.09019608,0.32941176,0.52941176>'
    chemical_color['Pt']='<0.81568627,0.81568627,0.87843137>'
    chemical_color['Au']='<1.00000000,0.81960784,0.13725490>'
    chemical_color['Hg']='<0.72156863,0.72156863,0.81568627>'
    chemical_color['Tl']='<0.65098039,0.32941176,0.30196078>'
    chemical_color['Pb']='<0.34117647,0.34901961,0.38039216>'
    chemical_color['Bi']='<0.61960784,0.30980392,0.70980392>'
    chemical_color['Po']='<0.67058824,0.36078431,0.00000000>'
    chemical_color['At']='<0.45882353,0.30980392,0.27058824>'
    chemical_color['Rn']='<0.25882353,0.50980392,0.58823529>'
    color=chemical_color[sym]
    return color
#------------------------------------------------------------------------------------------
