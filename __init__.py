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

    def fetch_components(self, comp=[0,1,2]):
        """ 
        This returns a graph of the components specified in num
        
        Args:
          comp: the components to return (default = top3)
            if list, returns those specific components, if number return range
        Return:
          A list of the components
        """
        mem = self.components().membership
        compcnt = {}
        for k in mem:
            if k not in compcnt:
                compcnt[k] = 0
            compcnt[k] = compcnt[k] + 1
        top = sorted(compcnt.items(), key=operator.itemgetter(1), reverse=True)

        if comp.__class__.__name__ not in "list":
            comp = range(0, comp)

        in_t = [top[c][0] for c in comp]
        return self.vs([i for i,m in enumerate(mem) if m in in_t]).subgraph()
