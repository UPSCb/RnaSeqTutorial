# -*- coding: utf-8 -*-
"""
Created on Tue May  5 17:43:19 2015

@author: bastianschiffthaler
"""

from sys import argv
import re

inf = argv[1]
ouf = argv[2]

OUF = open(ouf,"a")
ALL = []
with open(inf, "r") as INF:
    CODE=False
    BADCOMMANDS = [r"\\Biocpkg{(.+?)}",r"\\Rfunction{(.+?)}",
    r"\\Rclass{(.+?)}",r"\\Biocannopkg{(.+?)}",
    r"\\Rpackage{(.+?)}",r"\\bioccomment{(.+?)}",
    r"\\Rcode{(.+?)}",r"\\cite{(.+?)}"]
    REPLACEMENTS = {"\R{}" : "R",
                    "\Bioconductor{}" : "Bioconductor",
                    "\Dmel{}" : "\emph{Dmel}",
                    "\\begin{Exercise}" : "",
                    "\end{Exercise}" : "",
                    "\ie" : "\emph{i.e.}:",
                    "\eg" : "\emph{e.g.}:",
                    "\\begin{Solution}" : "",
                    "\end{Solution}" : "",
                    "\LaTeX{}":"LaTeX"} 
    
    for line in INF.readlines():
        if line.startswith("<<"):
            CODE=True
            ALL.append("\\begin{verbatim}")
            ALL.append("\n")
        elif line.startswith("@"):
            CODE=False
            ALL.append("\end{verbatim}")
        elif CODE:
            ALL.append("#" + line.replace("$","\$"))
        else:
            ALL.append(line)

ALL = "".join(ALL)

for i in REPLACEMENTS:
    ALL = ALL.replace(i,REPLACEMENTS[i])                
for f in BADCOMMANDS:
    ALL = re.sub(f,r'REPLACESQ' + r'\1' + "REPLACESQ",ALL,count=200)
OUF.write(ALL)

OUF.close()