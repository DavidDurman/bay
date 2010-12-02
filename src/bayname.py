#!/usr/bin/env python

# @project Learning Bayes network from data
# @author David Durman, 2008
# @description <p>This is a tool to rename nodes in a dot file
#               according to their real names.</p>

import sys

usage = """Usage: bayname.py namesfile dotfile
        namesfile       - file with node names in ascending order
                        - each line contains a node name
        dotfile         - Graphviz dot file generated from baybuild
"""

if __name__ == "__main__":
    try:
        namesFilename = sys.argv[1]
    except IndexError:
        sys.stderr.write(usage)
        sys.exit(1)

    try:
        dotFilename = sys.argv[2]
    except IndexError:
        sys.stderr.write(usage)
        sys.exit(1)

    namesFile = open(namesFilename, "r")
    dotFile = open(dotFilename, "r")

    names = []

    linenum = 0
    for line in namesFile.xreadlines():
        names.append(line.strip())
        linenum += 1

    namesFile.close()

    for line in dotFile.xreadlines():
        for i in xrange(linenum + 1):
            if line.find("\"node_" + str(i) + "\"") != -1:
                line = line.replace("\"node_" + str(i) + "\"", "\"" + names[i] + "\"")
        sys.stdout.write(line)
        
    dotFile.close()
