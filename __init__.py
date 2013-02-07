from igraph import *
from igraphex.drawing import *
from igraphex.drawing.graph import DefaultGraphDrawer2

class Graph(Graph):
    def __plot__(self, context, bbox, palette, *args, **kwds):
        """ 
        This overrides the original function but just includes DefaultGraphDrawer2
        """
        drawer_factory = kwds.get("drawer_factory", DefaultGraphDrawer2)
        if "drawer_factory" in kwds:
            del kwds["drawer_factory"]
        drawer = drawer_factory(context, bbox)
        drawer.draw(self, palette, *args, **kwds)

    def top_components(self, num=3, debug=False):
        """ 
        This allows us to quickly return n number of top components
        
        Args:
          num: number of components to return
        Return:
          A list of the components
        """
        mem = self.components().membership
        comp = {}
        for k in mem:
            if k not in comp:
                comp[k] = 0
            comp[k] = comp[k] + 1
        top = sorted(comp.items(), key=operator.itemgetter(1), reverse=True)[:num]

        c = []
        for t in top:
            c.append(self.vs(
                [i for i,m in enumerate(mem) if m==t[0]]).subgraph())
            if debug:
                print t, c[-1].vcount(), c[-1].ecount()
        return c

