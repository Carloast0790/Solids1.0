import os.path
bilatu_file='INPUT.txt'
#------------------------------------------------------------------------------------------
def read_block_of_bil(id):
    bilfile=open(bilatu_file,"r")
    chainchar='---'+id.upper()+'---'
    printer=0
    data_block=[]
    for line in bilfile:
         lin = line.lstrip()
         if lin.startswith(chainchar): printer=1+printer
         if printer == 1 and not lin.startswith(chainchar): data_block.append(line)
    bilfile.close()
    return data_block 
#------------------------------------------------------------------------------------------
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