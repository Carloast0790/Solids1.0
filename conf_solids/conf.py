import os
import os.path
from conf.bdata import *
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from conf.conf import host_name
    >>> HPC_NAME=host_name()
    >>> print("HPC_NAME = %s" %(HPC_NAME))
    """

def host_name():
    myhost = os.uname()[1]
    HPC_NAME = str(myhost)
    if HPC_NAME not in list_of_know_host: HPC_NAME='unknown'
    return HPC_NAME
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from conf.conf import queue_name
    >>> QUEUE_NAME=queue_name()
    >>> print("QUEUE_NAME = %s" %(QUEUE_NAME))
    """

def queue_name():
    HPC_NAME=host_name()
    queuename = default_queue.get(HPC_NAME)
    return queuename
#------------------------------------------------------------------------------------------
""" programg: 'gaussian', 'nwchem', 'vasp', 'mopac
    HOW TO USE:
    >>> from conf.conf import config_string
    >>> config_string('gaussian')
    """

def config_string(programg):
    HPC_NAME = host_name()
    print("HPC_NAME = %s" %(HPC_NAME))
    string_conf=long_string_conf[HPC_NAME, programg]
    return string_conf
#------------------------------------------------------------------------------------------
""" programg: 'gaussian', 'nwchem', 'vasp', 'mopac'
    HOW TO USE:
    >>> from conf.conf import config_file
    >>> config_file('gaussian')
    """

def config_file(programg):
    data=[]
    fileg=file_conf[programg]
    if not os.path.isfile(fileg):
        print("Creating configuration file :: %s, please check it" %(fileg))
        config=open(fileg,"w")
        string_conf=config_string(programg)
        config.write(string_conf)
        config.close()
    if os.path.isfile(fileg):
        fileconf=open(fileg,"r")
        for line in fileconf: data.append(line)
        fileconf.close()
    return data
#------------------------------------------------------------------------------------------
