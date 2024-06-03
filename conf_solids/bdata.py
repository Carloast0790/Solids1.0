#------------------------------------------------------------------------------------------
list_of_know_host=[
'kukulcan',
'dacb-hpc',
'login1',
'iMac-de-Filiberto.local',
'MacBook-Pro-de-Filiberto.local'
]
#------------------------------------------------------------------------------------------
default_queue = {
'kukulcan':'qgpu',
'dacb-hpc':'workq',
'login1':'CI',
'iMac-de-Filiberto.local':'JLO',
'MacBook-Pro-de-Filiberto.local':'JLO',
'unkonow':'queue_name'
}
#------------------------------------------------------------------------------------------
file_conf={}
file_conf['gaussian']='conf.gau'
file_conf['vasp']    ='conf.vasp'
file_conf['mopac']   ='conf.mop'
file_conf['nwchem']  ='conf.nw'
#------------------------------------------------------------------------------------------
long_string_conf={}
#------------------------------------------------------------------------------------------
###KUKULCAN
#------------------------------------------------------------------------------------------
long_string_conf['kukulcan', 'gaussian']= """export INT=i8
cd $PBS_O_WORKDIR

env > entorno-$PBS_JOBID.txt
export g16root=/LUSTRE/software/intel
source $g16root/g16/bsd/g16.profile
export PATH=$g16root/g16/:$PATH
export GAUSS_SCRDIR=/scratch/$USER/$PBS_JOBNAME-$PBS_JOBID
exe_gaussian=$g16root/g16/g16
mkdir -p $GAUSS_SCRDIR

RUNLINE

mv $GAUSS_SCRDIR $PBS_O_WORKDIR
rm -rf $PBS_JOBNAME-$PBS_JOBID
rm entorno-$PBS_JOBID.txt
"""

long_string_conf['kukulcan', 'vasp']= """
export VASP_PP=/LUSTRE/home/vasp/installdir/PBE/potpawPBE54
export PATH=/opt/intel/impi/4.1.1.036/intel64/bin/:$PATH
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64/:/opt/intel/lib/intel64/:$PATH
exe_mpirun=/opt/intel/impi/4.1.1.036/intel64/bin/mpirun
exe_vasp=/LUSTRE/software/intel/vasp/vasp_std-intel
"""

long_string_conf['kukulcan', 'mopac']="""
export LD_LIBRARY_PATH=/LUSTRE/home/fortiz/installdir/bin:$LD_LIBRARY_PATH
export MOPAC_LICENSE=/LUSTRE/home/fortiz/installdir/bin
exe_mopac=/LUSTRE/home/fortiz/installdir/bin/MOPAC2016.exe
"""
#------------------------------------------------------------------------------------------
###JUCHIMAN
#------------------------------------------------------------------------------------------
long_string_conf['dacb-hpc', 'gaussian']= """export INT=i8
export g16root=/home/apps/appg
source $g16root/g16/bsd/g16.profile
export PATH=$g16root/g16/:$PATH
export GAUSS_SCRDIR=$USER/$PBS_JOBNAME-$PBS_JOBID
exe_gaussian=$g16root/g16/g16
"""

long_string_conf['dacb-hpc', 'vasp']= """
export VASP_PP=/home/apps/appx/vasp544/PBE/potpawPBE54
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64/:/opt/intel/lib/intel64/:$PATH
export PATH=/opt/intel/impi/5.1.2.150/intel64/bin/:$PATH
exe_mpirun=/opt/intel/impi/5.1.2.150/intel64/bin/mpirun
exe_vasp=/home/apps/appx/vasp544/vasp_gam
"""

long_string_conf['dacb-hpc', 'mopac']="""
export MOPAC_LICENSE=/home/apps/appgnu/MOPAC/centos/
exe_mopac=/home/apps/appgnu/MOPAC/centos/MOPAC2016.exe
"""

long_string_conf['dacb-hpc', 'nwchem']= """
export PATH=/home/apps/nwchem-6.6/bin/LINUX64/:$PATH
export PATH=/opt/intel/impi/5.1.2.150/intel64/bin/:$PATH
source /home/apps/nwchem-6.6/nwchem.vars
NWCHEM_BASIS_LIBRARY=/home/apps/nwchem-6.6/src/basis/libraries
exe_nwchem=/home/apps/nwchem-6.6/bin/LINUX64/nwchem
exe_mpirun=/opt/intel/impi/5.1.2.150/intel64/bin/mpirun
"""
#------------------------------------------------------------------------------------------
###Xiuhcoatl
#------------------------------------------------------------------------------------------
long_string_conf['login1', 'gaussian']= """
cd $PBS_O_WORKDIR

env > entorno-$PBS_JOBID.txt
export g09root=/usr/local/g09ID.01
source $g09root/g09/bsd/g09.profile
export GAUSS_SCRDIR=/scratch/users/$USER/$PBS_JOBNAME-$PBS_JOBID
exe_gaussian=$g09root/g09/g09
mkdir -p $GAUSS_SCRDIR

RUNLINE

mv $GAUSS_SCRDIR $PBS_O_WORKDIR
rm -rf $PBS_JOBNAME-$PBS_JOBID
rm entorno-$PBS_JOBID.txt
rm $PBS_JOBNAME.pbs
"""
#------------------------------------------------------------------------------------------
###iMac-de-Filiberto.local
#------------------------------------------------------------------------------------------
long_string_conf['iMac-de-Filiberto.local', 'mopac']="""
export MOPAC_LICENSE=/Users/fortiz/installdir/bin
exe_mopac=/Users/fortiz/installdir/bin/MOPAC2016.exe
"""

long_string_conf['iMac-de-Filiberto.local', 'vasp']= """
export VASP_PP=/Users/fortiz/installdir/POT.VASP/potPAW.54_PBE
export PATH=/opt/intel/impi/5.1.2.150/intel64/bin/:$PATH
exe_mpirun=/opt/intel/impi/5.1.2.150/intel64/bin/mpirun
exe_vasp=/home/apps/appx/vasp/vasp.5.4.1_intel/bin/vasp_gam
"""
#------------------------------------------------------------------------------------------
###MacBook-Pro-de-Filiberto.local
#------------------------------------------------------------------------------------------
long_string_conf['MacBook-Pro-de-Filiberto.local', 'mopac']="""
export MOPAC_LICENSE=/Users/fortiz/installdir/bin
exe_mopac=/Users/fortiz/installdir/bin/MOPAC2016.exe
"""
#------------------------------------------------------------------------------------------
###UNKNOWN
#------------------------------------------------------------------------------------------
long_string_conf['unkonow', 'gaussian']= """
export g09root=/opt/intel
source $g09root/g09/bsd/g09.profile
export PATH=$PATH:/opt/intel/g09/g09
export GAUSS_SCRDIR=/scratch/$USER/$PBS_JOBNAME-$PBS_JOBID
"""

long_string_conf['unkonow', 'mopac']= """
export PATH=$HOME/installdir/bin/:$PATH
export MOPAC_LICENSE=$HOME/installdir/bin
exe_mopac=$HOME/installdir/bin/MOPAC2016.exe
"""
#------------------------------------------------------------------------------------------
