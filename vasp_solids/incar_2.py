long_string = """general:
 SYSTEM  = Generic !Identify what to do with this specific input file
 PREC = Normal     !Standard precision 
 ISTART  = 0       !Startjob: 0-new 1-cont 2-samecut
 ENMAX = 400       !Sutoff should be set manually
 ICHARG = 2        !initial charge density: 1-file 2-atom 10-cons 11-DOS
 ISMEAR = 0        !k-mesh integration: 0 Gaussian smearing
 SIGMA = 0.02      !Insulators/semiconductors=0.1  metals=0.05
 IBRION = 2        !Use DIIS algorithm to converge
 LREAL = .FALSE.   !Projection done in reciprocal space
 ADDGRID = .TRUE   !Determines whether an additional support grid is used for the evaluation of the augmentation
 charges
 NFREE = 2         !Determines how many displacements are used for each direction and ion
 NSW = 150         !Number of ionic steps
 NELM = 60         !maximum number of electronic SC
 NELMDL = 6        !6 steps with Davidson; anothe with the default
 EDIFFG = -0.02    !forces smaller 0.02 A/eV
 ALGO=VeryFast     !Selects a faily robust mixture of the Davidson and RMM-DIIS algorithms
 LSCALAPACK = .FALSE.
spin:
 ISPIN   = 1       !Spin polarized calculation: 2-yes 1-no
 LCHARG = FALSE    !Avoid, that CHG and CHGCAR is written out
 LWAVE = FALSE     !Avoid, that WAVECAR is written out
"""

exfile = open("incar_2", "w")
exfile.write(long_string)
exfile.close()
