long_string = """Gamma-point only
 0
Monkhorst Pack
 1 1 1
 0 0 0
"""

exfile = open("kpoints", "w")
exfile.write(long_string)
exfile.close()
