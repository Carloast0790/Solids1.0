long_string = """Automatic kpoint scheme
0
Monkhorst
4 4 4
"""

exfile = open("kpoints", "w")
exfile.write(long_string)
exfile.close()
