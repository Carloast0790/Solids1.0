import os.path
import sys
from inout_solids.messag import write_welcome, write_menu
from utils_solids.atomic import check_atomic_sym, is_number, get_atomic_mass
#------------------------------------------------------------------------------------------
def get_bilfile():
    bilfile='INPUT.txt'
    if not os.path.isfile(bilfile):
        args=sys.argv[1:]
        if len(args) < 1:
            write_welcome()
            write_menu()
            print("INPUT.txt was not found!")
            print("   ")
            print("To set up your calculation, build an INPUT.txt file,")
            print("similar to what you want to do, and edit it. Try:")
            print("   ")
            print("solids --example_SM_vasp_TiO2")
            print("solids --example_SM_gulp_TiO2")
            print("solids --example_SM_gulp_C8")
            print("solids --example_SM_gulp_Si8")
            print("solids --example_GA_vasp_TiO2")
            print("solids --example_GA_gulp_TiO2")
            print("solids --example_GA_gulp_MgAl2O4")
            print("solids --example_GA_gulp_MgSiO3")
            print("solids --example_GA_gulp_SrTiO3")
            print("   ")
            sys.exit()
        #====================================================================
        ##VASP
        #====================================================================
        if sys.argv[1] == "--example_SM_vasp_TiO2":
            print("Built INPUT.txt file for TiO2 search using VASP and the Modified Stochastic Method")
            import vasp_solids.exampleKV
            import vasp_solids.incar_1
            import vasp_solids.incar_2
            import vasp_solids.kpoints
            sys.exit()
        if sys.argv[1] == "--example_GA_vasp_TiO2":
            print("Built INPUT.txt file for TiO2 search using VASP and the Modified Genetic Algorithm")
            import vasp_solids.exampleGAV
            import vasp_solids.incar_1
            import vasp_solids.incar_2
            import vasp_solids.kpoints
            sys.exit()
        #====================================================================
        ##GULP
        #====================================================================
        if sys.argv[1] == "--example_SM_gulp_TiO2":
            print("Built INPUT.txt file for TiO2 search using GULP and the Modified Stochastic Method")
            import gulp_solids.example_SM_gulp_TiO2

        if sys.argv[1] == "--example_SM_gulp_C8":
            print("Built INPUT.txt file for C8 search using GULP and the Modified Stochastic Method")
            import gulp_solids.example_SM_gulp_C8

        if sys.argv[1] == "--example_SM_gulp_Si8":
            print("Built INPUT.txt file for Si8 search using GULP and the Modified Stochastic Method")
            import gulp_solids.example_SM_gulp_Si8

        if sys.argv[1] == "--example_GA_gulp_TiO2":
            print("Built INPUT.txt file for TiO2 search using GULP and the Modified Genetic Algorithm")
            import gulp_solids.example_GA_gulp_TiO2

        if sys.argv[1] == "--example_GA_gulp_MgAl2O4":
            print("Built INPUT.txt file for MgAl2O4 search using GULP and the Modified Genetic Algorithm")
            import gulp_solids.example_GA_gulp_MgAl2O4

        if sys.argv[1] == "--example_GA_gulp_MgSiO3":
            print("Built INPUT.txt file for MgSiO3 search using GULP and the Modified Genetic Algorithm")
            import gulp_solids.example_GA_gulp_MgSiO3

        if sys.argv[1] == "--example_GA_gulp_SrTiO3":
            print("Built INPUT.txt file for SrTiO3 search using GULP and the Modified Genetic Algorithm")
            import gulp_solids.example_GA_gulp_SrTiO3

        print("Edit the INPUT.txt file and re-run solids as:")
        print("nohup solids > log &")
        sys.exit()
    else:
        return bilfile
#------------------------------------------------------------------------------------------
def read_block_of_bil(id):
    bilatu_file=get_bilfile()
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
def read_var_composition(id):
    bilatu_file=get_bilfile()
    bilfile=open(bilatu_file,"r")
    chainchar='---'+id.upper()+'---'
    printer=0
    data_block=[]
    for line in bilfile:
         lin = line.lstrip()
         if lin.startswith(chainchar): printer=1+printer
         if printer == 1 and not lin.startswith(chainchar):
             data_block.append(line.split())
    bilfile.close()
    nat, nvc, varcomp=len(data_block), len(data_block[0]), []
    for jj  in range(1,nvc):
        icomp=[]
        for ii in range(nat):
            atom=data_block[ii][0]
            nats=int(data_block[ii][jj])
            icomp.append([atom,nats])
        varcomp.append(icomp)
    return varcomp
#------------------------------------------------------------------------------------------
def clustername(composition):
    chainname=' '.join([item[0]+str(item[1]) for item in composition])
    return chainname
#------------------------------------------------------------------------------------------
def get_inatoms(composition):
    inatoms=[]
    for xxx in composition:
        for jj in range(xxx[1]):
            inatoms.append(xxx[0])
    return inatoms
#------------------------------------------------------------------------------------------
#QUEREMOS ERRADICAR ESTE:
def get_stchiom(opt=0):
    atoms=[]
    nofatoms=[]
    bilatu_file=get_bilfile()
    bilfile=open(bilatu_file,"r")
    for line in bilfile:
        line=line.strip(' \t\n\r')
        if len(line.strip()) != 0 :
            li = line.lstrip()
            if not li.startswith("#"):
                readline=line.split()
                if len(readline) == 2:
                    sym=readline[0].strip(' \t\n\r')
                    nato=readline[1].strip(' \t\n\r')
                    if check_atomic_sym(sym) and is_number(nato):
                        atoms.append(sym)
                        nofatoms.append(int(nato))
    bilfile.close()
    listm=[get_atomic_mass(x) for x in atoms]
    s = list(zip(atoms, nofatoms, listm))
    t = sorted(s, key=lambda x: float(x[2]), reverse=True)

    allatoms=([(x[0], x[1]) for x in t])
    atoms=([x[0] for x in t])       
    nofatoms=([x[1] for x in t])
    clustername=' '.join([item[0]+str(item[1]) for item in allatoms])
    inatoms=[]
    for iii in allatoms:
        for jjj in range(iii[1]):   
            inatoms.append(iii[0])
    return {
	0: allatoms[:],
        1: atoms,
        2: nofatoms,
        3: clustername,
        4: inatoms
    }[opt]
#------------------------------------------------------------------------------------------
