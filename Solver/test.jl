using Graphs
using LightGraphs

const _GR = Graphs
const _LG = LightGraphs

g = _LG.DiGraph(4)

_LG.add_edge!(g, 1, 2)
_LG.add_edge!(g, 2, 3)
_LG.add_edge!(g, 3, 4)
_LG.add_edge!(g, 4, 1)
_LG.add_edge!(g, 1, 3)

_LG.outneighbors(g,1)