solids_short_string="""# Solids 1.0
"""

welcome_solids = """# Solids 1.0
#
# Authors:
# (1) Carlos López-Castro, (2) Filiberto Ortiz-Chi, and (1) Gabriel Merino
#
# (1) Departamento de Fisica Aplicada,
#     Centro de Investigacion y de Estudios Avanzados, Unidad Merida,
#     Km 6 Antigua Carretera a Progreso, A.P. 73, Cordemex, 97310 Merida, Yucatan, Mexico
# (2) CONACYT-Universidad Juarez Autonoma de Tabasco,
#     Centro de Investigacion de Ciencia y Tecnologia Aplicada de Tabasco,
#     Cunduacan 86690, Tabasco, Mexico.
"""


menu_long_string = """
#----------------------------------------------------------------------------
# SYSTEM  | METHOD                     |    CODE   |  opt  |  NEEDS
#----------------------------------------------------------------------------
#Crystal  | Modified Stochastic Method | VASP/GULP | kick  | (COMPOSITION)
#Crystal  | Modified Genetic Algorithm | VASP/GULP |  GA   | (COMPOSITION)
#----------------------------------------------------------------------------
"""

def write_welcome(logfile=''):
    if logfile=='':
        print(welcome_solids)
    else:
        print(welcome_solids, file=logfile)

def write_menu():
    print(menu_long_string)

def head_gega():
    print(welcome_glomos)

def head_growpal():
    print(welcome_growpal)

def head_solids():
    print(welcome_solids)