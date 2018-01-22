def addScaffInfo(filename, f, color, name, idbyname, cntlib, cntnode, cntedge):
    with open(filename) as g:
        f.write("scaffoldlibs.push(new ScaffoldEdgeLib(" + str(cntlib) + ", '" + str(color) + "', '" + str(name) + "', 'SCAFF'));\n")
        cntlib += 1

        scafnum = 0

        for line in g:
            tokens = line.split(" ")
            if (tokens[len(tokens) - 1] == '\n'):
                tokens.pop()

            f.write("scaffoldlibs["+ str(cntlib - 1) +"].scaffolds.push(new Scaffold('" + tokens[0][1:] + "'));\n")
            f.write("scaffoldlibs["+ str(cntlib - 1) +"].scaffolds.push(new Scaffold('" + tokens[0][1:] + "-rev'));\n")

            nodeslist = []
            for i in range(1, len(tokens), 3):
                nm = tokens[i][1:]
                if (tokens[i + 2][0] == '+'):
                    nodeslist.append(idbyname[nm])
                else:
                    nodeslist.append(idbyname[nm]^1)


            for i in range(1, len(nodeslist)):
                f.write("scaffoldedges.push(new ScaffoldEdge(" + str(cntedge) + ", "+ str(nodeslist[i - 1]) +
                        ", " + str(nodeslist[i]) + ", " + str(cntlib - 1) + ", 1));\n")
                f.write("scaffoldedges["+str(cntedge)+"].name='"+ tokens[0][1:] + "';\n")
                f.write("scaffoldlibs["+ str(cntlib - 1) +"].scaffolds["+str(scafnum) +"].edges.push(scaffoldedges["+str(cntedge)+"]);\n")
                cntedge += 1

            for i in range(len(nodeslist)-2, -1, -1):
                f.write("scaffoldedges.push(new ScaffoldEdge(" + str(cntedge) + ", "+ str(nodeslist[i + 1]^1) +
                        ", " + str(nodeslist[i]^1) + ", " + str(cntlib - 1) + ", 1));\n")
                f.write("scaffoldedges["+str(cntedge)+"].name='"+ tokens[0][1:] + "-rev';\n")
                f.write("scaffoldlibs["+ str(cntlib - 1) +"].scaffolds["+str(scafnum+1) +"].edges.push(scaffoldedges["+str(cntedge)+"]);\n")
                cntedge += 1

            scafnum += 2

    return cntlib, cntnode, cntedge




g = open("data/graph.gr", "r")
f = open("genescaffoldgraphb.js", "w")

idbyname = dict()

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
    idbyname[nodesinfo[2]] = i
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

libcnt, nodecnt, edgescnt = addScaffInfo("data/out.info", f, "#ffff00", "IGAR", idbyname, libcnt, nodecnt, edgescnt)
libcnt, nodecnt, edgescnt = addScaffInfo("data/rascaf.info", f, "#ff00ff", "rascaf", idbyname, libcnt, nodecnt, edgescnt)

f.write("var scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);")

f.close()
g.close()

af = open("genealignmentb.js", "w")

chrlist = []
chralig = []

curid = -2
lastname = '-'

with open("data/out.coords") as cf:
    for line in cf:
        info = line.split("\t")
        info[10] = info[10][:-1]
        vid = idbyname[info[10]]
        chrname = info[9]
        lenf = int(info[7])
        if (chrname != lastname):
            curid += 2
            chrlist.append("new Chromosome(" + str(curid) + ", '" + chrname + "', " + str(lenf) + ")")
            chrlist.append("new Chromosome(" + str(curid + 1) + ", '" + chrname + "-rev', " + str(lenf) + ")")
            chralig.append([])
            chralig.append([])
            lastname = chrname

        lq = int(info[2])
        rq = int(info[3])
        l = int(info[0])
        r = int(info[1])

        if (lq > rq):
            vid ^= 1

        chralig[curid].append("new Alignment(" + str(l) + ", " + str(r) + ", " + str(curid) + ", " + str(vid) + ")")
        chralig[curid + 1].append("new Alignment(" + str(lenf - r) + ", " + str(lenf - l) + ", " + str(curid + 1) + ", " + str(vid^1) + ")")

af.write("var chromosomes = [")
for i in range(len(chrlist)):
    af.write(chrlist[i])
    if (i != len(chrlist) - 1):
        af.write(", ")
af.write("];\n")

for i in range(len(chrlist)):
    af.write("chromosomes[" + str(i) + "].alignments = " + "[")
    for j in range(len(chralig[i])):
        af.write(chralig[i][j])
        if (j != len(chralig[i]) - 1):
            af.write(", ")
    af.write("];\n")

af.close()