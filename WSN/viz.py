import numpy as np
import toyplot
import toyplot.png
import toyplot.pdf
import sys

graph = open("graph.txt","r")

line1 = graph.readline()
line2 = graph.readline()
line3 = graph.readline()

edges = line2.split(" ")
edges.pop()
edges = [int(x) for x in edges]

offset = line3.split(" ")
offset.pop()
offset = [int(x) for x in offset]

graph.close()

cords = open("Cord.txt","r")

xy = []
for i in range(0,len(offset)):
    tline = cords.readline()
    tline = tline.split(" ")
    tline.pop()
    tline = [float(x) for x in tline]
    xy.append(tline)

cords.close()
xy = np.array(xy)

extra_nodes=[]
tedges = []
for i in range(0,len(offset)-1): 
    start = offset[i]
    end = offset[i+1]
    if start == end:
        extra_nodes.append(i)
    while start < end:
        l = [i,edges[start]]
        tedges.append(l)
        start+=1


vcoordinates = np.ma.masked_all((len(offset), 2)) 
for i in range(0,len(offset)):
    vcoordinates[i]= ( xy[i][0], xy[i][1])


colormap = toyplot.color.LinearMap(toyplot.color.Palette(["white", "yellow", "red"]))
vstyle = {"stroke":toyplot.color.black}

tedges = np.array(tedges)
layout = toyplot.layout.FruchtermanReingold()
c,y,z = toyplot.graph(tedges,extra_nodes,vcoordinates=vcoordinates,width=600)#,vcolor=colormap, vsize=5, vstyle=vstyle)
# toyplot.png.render(c, "figure1.png")
toyplot.pdf.render(c, "figure1.pdf")

c1,y1,z1 = toyplot.graph(tedges,width=600)#,vcolor=colormap, vsize=5, vstyle=vstyle)
toyplot.pdf.render(c1, "figure2.pdf")