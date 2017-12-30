g = open("graph.gr", "r")
f = open("genescaffoldgraph.js", "w")

libstr = "var scaffoldlibs = ["

libcnt = int(g.readline())
for i in range(libcnt):
    libsinfo = g.readline().split(" ")
    libsinfo[4] = libsinfo[4][:-1]
    if (i != 0):
        libstr += ', '
    libstr += "new ScaffoldEdgeLib(" + libsinfo[1] + ", '" + libsinfo[2] + "', '" + libsinfo[3] + "', '" + libsinfo[4] + "')"

libstr += "];"

nodestr = "var scaffoldnodes = ["

nodecnt = int(g.readline())
for i in range(nodecnt):
    nodesinfo = (g.readline()).split(" ")
    nodesinfo[3] = nodesinfo[3][:-1]
    if (i != 0):
        nodestr += ', '
    nodestr += "new ScaffoldNode(" + nodesinfo[1] + ", '" + nodesinfo[2] + "', " + nodesinfo[3] + ")"

nodestr += "];"

edgestr = "var scaffoldedges = ["

edgescnt = int(g.readline())
for i in range(edgescnt):
    edgesinfo = (g.readline()).split(" ")
    if (i != 0):
        edgestr += ', '
    edgestr += "new ScaffoldEdge(" + edgesinfo[1] + ", " + edgesinfo[2] + ", " + edgesinfo[3] + ", " + edgesinfo[4] + ", " + edgesinfo[5] + ")"

edgestr += "];"

f.write(libstr + "\n")
f.write(nodestr + "\n")
f.write(edgestr + "\n")
f.write("var scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);")

f.close()
g.close()

