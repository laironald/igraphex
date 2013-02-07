from igraph import *
from igraphex.drawing import *
from igraphex.drawing.graph import DefaultGraphDrawer2

class Graph(Graph):
    def __plot__(self, context, bbox, palette, *args, **kwds):
        drawer_factory = kwds.get("drawer_factory", DefaultGraphDrawer2)
        if "drawer_factory" in kwds:
            del kwds["drawer_factory"]
        drawer = drawer_factory(context, bbox)
        drawer.draw(self, palette, *args, **kwds)

