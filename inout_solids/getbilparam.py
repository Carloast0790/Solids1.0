import os.path
bilatu_file='INPUT.txt'
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.getbilparam import get_a_float
    >>> b=get_a_float('convpercent',80)
    """

def get_a_float(strchain, defaultvalue):
    finalvalue=float(defaultvalue)
    if os.path.isfile(bilatu_file):
        bilfile=open(bilatu_file,"r")
        for line in bilfile:
            line=line.strip(' \t\n\r')
            if len(line.strip()) != 0 :
                li = line.lstrip()
                if not li.startswith("#"):
                    readline=line.split()
                    if len(readline) == 2:
                        data0=readline[0].strip('\t\n\r') 
                        data1=readline[1].strip('\t\n\r')
                        if data0.lower() == str(strchain): finalvalue=float(data1)
        bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.getbilparam import get_a_int
    >>> b=get_a_int('njobs',10)
    """

def get_a_int(strchain, defaultvalue):
    finalvalue=int(defaultvalue)
    if os.path.isfile(bilatu_file):
        bilfile=open(bilatu_file,"r")
        for line in bilfile:
            line=line.strip(' \t\n\r')
            if len(line.strip()) != 0 :
                li = line.lstrip()
                if not li.startswith("#"):
                    readline=line.split()
                    if len(readline) == 2:
                        data0=readline[0].strip('\t\n\r') 
                        data1=readline[1].strip('\t\n\r')
                        if data0.lower() == str(strchain): finalvalue=int(data1)
        bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.getbilparam import get_a_str
    >>> b=get_a_str('walltime','02:00:00')
    """

def get_a_str(strchain, defaultvalue):
    finalvalue=str(defaultvalue)
    if os.path.isfile(bilatu_file):
        bilfile=open(bilatu_file,"r")
        for line in bilfile:
            line=line.strip(' \t\n\r')
            if len(line.strip()) != 0 :
                li = line.lstrip()
                if not li.startswith("#"):
                    readline=line.split()
                    if len(readline) == 2:
                        data0=readline[0].strip('\t\n\r') 
                        data1=readline[1].strip('\t\n\r')
                        if data0.lower() == str(strchain): finalvalue=str(data1)
        bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.getbilparam import get_float_list
    >>> b=get_float_list('percent_of_convergence',[90.0, 90.0])
    """

def get_float_list(strchain, defaultvalue):
    finalvalue=[float(item) for item in defaultvalue]
    if os.path.isfile(bilatu_file):
        bilfile=open(bilatu_file,"r")
        for line in bilfile:
            line=line.strip(' \t\n\r')
            if len(line.strip()) != 0 :
                li = line.lstrip()
                if not li.startswith("#"):
                    readline=line.split()
                    data0=readline[0].strip('\t\n\r') 
                    if data0.lower() == str(strchain):
                        del readline[0]
                        data1=[float(item) for item in readline]
                        finalvalue=data1
        bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.getbilparam import get_int_list
    >>> b=get_int_list('num_of_attemps_to_conv',[2, 2])
    """

def get_int_list(strchain, defaultvalue):
    finalvalue=[int(item) for item in defaultvalue]
    if os.path.isfile(bilatu_file):
        bilfile=open(bilatu_file,"r")
        for line in bilfile:
            line=line.strip(' \t\n\r')
            if len(line.strip()) != 0 :
                li = line.lstrip()
                if not li.startswith("#"):
                    readline=line.split()
                    data0=readline[0].strip('\t\n\r') 
                    if data0.lower() == str(strchain):
                        del readline[0]
                        data1=[int(item) for item in readline]
                        finalvalue=data1
        bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.getbilparam import get_str_list
    >>> b=get_str_list('num_of_attemps_to_conv',['04:00:00', '12:00:00'])
    """

def get_str_list(strchain, defaultvalue):
    finalvalue=[str(item) for item in defaultvalue]
    if os.path.isfile(bilatu_file):
        bilfile=open(bilatu_file,"r")
        for line in bilfile:
            line=line.strip(' \t\n\r')
            if len(line.strip()) != 0 :
                li = line.lstrip()
                if not li.startswith("#"):
                    readline=line.split()
                    data0=readline[0].strip('\t\n\r') 
                    if data0.lower() == str(strchain):
                        del readline[0]
                        data1=[str(item) for item in readline]
                        finalvalue=data1
        bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
