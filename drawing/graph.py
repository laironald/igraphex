from igraph.drawing.graph import DefaultGraphDrawer

#-------------------------------------------
from collections import defaultdict
from itertools import izip
from math import atan2, cos, pi, sin, tan
from warnings import warn

from igraph._igraph import convex_hull, VertexSeq
from igraph.configuration import Configuration
from igraph.drawing.baseclasses import AbstractDrawer, AbstractCairoDrawer, \
                                       AbstractXMLRPCDrawer
from igraph.drawing.colors import color_to_html_format, color_name_to_rgb
from igraph.drawing.edge import ArrowEdgeDrawer
from igraph.drawing.text import TextDrawer
from igraph.drawing.metamagic import AttributeCollectorBase, \
                                     AttributeSpecification
from igraph.drawing.shapes import ShapeDrawerDirectory, PolygonDrawer
from igraph.drawing.utils import Point
from igraph.layout import Layout

__all__ = ["DefaultGraphDrawer", "UbiGraphDrawer", "CytoscapeGraphDrawer"]
__license__ = "GPL"

try:
    import cairo
except ImportError:
    # No cairo support is installed. Create a fake module
    # pylint: disable-msg=C0103
    from igraph.drawing.utils import FakeModule
    cairo = FakeModule()
#-------------------------------------------

class DefaultGraphDrawer2(DefaultGraphDrawer):
    def draw(self, graph, palette, *args, **kwds):
        # Some abbreviations for sake of simplicity
        directed = graph.is_directed()
        context = self.context

        # Calculate/get the layout of the graph
        layout = self.ensure_layout(kwds.get("layout", None), graph)

        # Determine the size of the margin on each side
        margin = kwds.get("margin", 0)
        try:
            margin = list(margin)
        except TypeError:
            margin = [margin]
        while len(margin)<4:
            margin.extend(margin)
        # margin = [x + 20. for x in margin[:4]]

        # Contract the drawing area by the margin and fit the layout
        bbox = self.bbox.contract(margin)
        layout.fit_into(bbox, keep_aspect_ratio=False)

        # Decide whether we need to calculate the curvature of edges
        # automatically -- and calculate them if needed.
        autocurve = kwds.get("autocurve", None)
        if autocurve or (autocurve is None and \
                "edge_curved" not in kwds and "curved" not in graph.edge_attributes() \
                and graph.ecount() < 10000):
            from igraph import autocurve
            default = kwds.get("edge_curved", 0)
            if default is True:
                default = 0.5
            default = float(default)
            kwds["edge_curved"] = autocurve(graph, attribute=None, default=default)

        # Construct the visual vertex/edge builders
        class VisualVertexBuilder(AttributeCollectorBase):
            """Collects some visual properties of a vertex for drawing"""
            _kwds_prefix = "vertex_"
            color = ("green", palette.get)
            label = None
            label_angle = -pi/2
            label_dist  = 0.0
            label_color = ("black", palette.get)
            label_size  = 14.0
            position = dict(func=layout.__getitem__)
            shape = ("circle", ShapeDrawerDirectory.resolve_default)
            size  = 20.0
            # RON EDIT ---------------------
            edge_color = ("#444", palette.get)
            edge_size = 2.0
            # ------------------------------

        class VisualEdgeBuilder(AttributeCollectorBase):
            """Collects some visual properties of an edge for drawing"""
            _kwds_prefix = "edge_"
            arrow_size  = 1.0
            arrow_width = 1.0
            color       = ("#444", palette.get)
            curved      = (0.0, ArrowEdgeDrawer._curvature_to_float)
            width       = 1.0

        vertex_builder = VisualVertexBuilder(graph.vs, kwds)
        edge_builder = VisualEdgeBuilder(graph.es, kwds)

        # Draw the highlighted groups (if any)
        if "mark_groups" in kwds:
            mark_groups = kwds["mark_groups"]

            # Figure out what to do with mark_groups in order to be able to
            # iterate over it and get memberlist-color pairs
            if isinstance(mark_groups, dict):
                group_iter = mark_groups.iteritems()
            elif hasattr(mark_groups, "__iter__"):
                # Lists, tuples, iterators etc
                group_iter = iter(mark_groups)
            else:
                # False
                group_iter = {}.iteritems()

            # We will need a polygon drawer to draw the convex hulls
            polygon_drawer = PolygonDrawer(context, bbox)

            # Iterate over color-memberlist pairs
            for group, color_id in group_iter:
                if not group or color_id is None:
                    continue

                color = palette.get(color_id)

                if isinstance(group, VertexSeq):
                    group = [vertex.index for vertex in group]
                if not hasattr(group, "__iter__"):
                    raise TypeError("group membership list must be iterable")

                # Get the vertex indices that constitute the convex hull
                hull = [group[i] for i in convex_hull([layout[idx] for idx in group])]

                # Calculate the preferred rounding radius for the corners
                corner_radius = 1.25 * max(vertex_builder[idx].size for idx in hull)

                # Construct the polygon
                polygon = [layout[idx] for idx in hull]

                if len(polygon) == 2:
                    # Expand the polygon (which is a flat line otherwise)
                    a, b = Point(*polygon[0]), Point(*polygon[1])
                    c = corner_radius * (a-b).normalized()
                    n = Point(-c[1], c[0])
                    polygon = [a + n, b + n, b - c, b - n, a - n, a + c]
                else:
                    # Expand the polygon around its center of mass
                    center = Point(*[sum(coords) / float(len(coords))
                                      for coords in zip(*polygon)])
                    polygon = [Point(*point).towards(center, -corner_radius)
                               for point in polygon]

                # Draw the hull
                context.set_source_rgba(color[0], color[1], color[2],
                                        color[3]*0.25)
                polygon_drawer.draw_path(polygon, corner_radius=corner_radius)
                context.fill_preserve()
                context.set_source_rgba(*color)
                context.stroke()

        # Draw the edges
        edge_drawer = self.edge_drawer_factory(context)
        if directed:
            drawer_method = edge_drawer.draw_directed_edge
        else:
            drawer_method = edge_drawer.draw_undirected_edge
        for edge, visual_edge in izip(graph.es, edge_builder):
            src, dest = edge.tuple
            src_vertex, dest_vertex = vertex_builder[src], vertex_builder[dest]
            drawer_method(visual_edge, src_vertex, dest_vertex)

        # Calculate the desired vertex order
        if "vertex_order" in kwds:
            # Vertex order specified explicitly
            vertex_order = kwds["vertex_order"]
        elif kwds.get("vertex_order_by") is not None:
            # Vertex order by another attribute
            vertex_order_by = kwds["vertex_order_by"]
            if isinstance(vertex_order_by, tuple):
                vertex_order_by, reverse = vertex_order_by
                if isinstance(reverse, basestring) and reverse.lower().startswith("asc"):
                    reverse = False
                else:
                    reverse = bool(reversed)
            else:
                reverse = False
            attrs = graph.vs[vertex_order_by]
            vertex_order = sorted(range(graph.vcount()), key=attrs.__getitem__,
                    reverse=reverse)
            del attrs
        else:
            # Default vertex order
            vertex_order = None

        if vertex_order is None:
            # Default vertex order
            vertex_coord_iter = izip(vertex_builder, layout)
        else:
            # Specified vertex order
            vertex_coord_iter = ((vertex_builder[i], layout[i])
                    for i in vertex_order)

        # RON EDIT ---------------------
        # Draw the vertices
        context.set_line_width(1)
        for vertex, coords in vertex_coord_iter:
            vertex.shape.draw_path(context, coords[0], coords[1], vertex.size)
            context.set_source_rgba(*vertex.color)
            context.fill_preserve()
            #print vertex.edge_color
            context.set_source_rgba(*vertex.edge_color)
            context.set_line_width(vertex.edge_size)
            context.stroke()
        # ------------------------------

        # Draw the vertex labels
        context.select_font_face("sans-serif", cairo.FONT_SLANT_NORMAL, \
            cairo.FONT_WEIGHT_NORMAL)

        wrap = kwds.get("wrap_labels", None)
        if wrap is None:
            wrap = Configuration.instance()["plotting.wrap_labels"]
        else:
            wrap = bool(wrap)

        if vertex_order is None:
            # Default vertex order
            vertex_coord_iter = izip(vertex_builder, layout)
        else:
            # Specified vertex order
            vertex_coord_iter = ((vertex_builder[i], layout[i])
                    for i in vertex_order)

        label_drawer = self.label_drawer_factory(context)
        for vertex, coords in vertex_coord_iter:
            if vertex.label is None:
                continue

            context.set_font_size(vertex.label_size)
            context.set_source_rgba(*vertex.label_color)
            label_drawer.text = vertex.label

            if vertex.label_dist:
                # Label is displaced from the center of the vertex.
                _, yb, w, h, _, _ = label_drawer.text_extents()
                w, h = w/2.0, h/2.0
                radius = vertex.label_dist * vertex.size / 2.
                # First we find the reference point that is at distance `radius'
                # from the vertex in the direction given by `label_angle'.
                # Then we place the label in a way that the line connecting the
                # center of the bounding box of the label with the center of the
                # vertex goes through the reference point and the reference
                # point lies exactly on the bounding box of the vertex.
                alpha = vertex.label_angle % (2*pi)
                cx = coords[0] + radius * cos(alpha)
                cy = coords[1] - radius * sin(alpha)
                # Now we have the reference point. We have to decide which side
                # of the label box will intersect with the line that connects
                # the center of the label with the center of the vertex.
                if w > 0:
                    beta = atan2(h, w) % (2*pi)
                else:
                    beta = pi/2.
                gamma = pi - beta
                if alpha > 2*pi-beta or alpha <= beta:
                    # Intersection at left edge of label
                    cx += w
                    cy -= tan(alpha) * w
                elif alpha > beta and alpha <= gamma:
                    # Intersection at bottom edge of label
                    try:
                        cx += h / tan(alpha)
                    except:
                        pass    # tan(alpha) == inf
                    cy -= h
                elif alpha > gamma and alpha <= gamma + 2*beta:
                    # Intersection at right edge of label
                    cx -= w
                    cy += tan(alpha) * w
                else:
                    # Intersection at top edge of label
                    try:
                        cx -= h / tan(alpha)
                    except:
                        pass    # tan(alpha) == inf
                    cy += h
                # Draw the label
                label_drawer.draw_at(cx-w, cy-h-yb, wrap=wrap)
            else:
                # Label is exactly in the center of the vertex
                cx, cy = coords
                half_size = vertex.size / 2.
                label_drawer.bbox = (cx - half_size, cy - half_size,
                                     cx + half_size, cy + half_size)
                label_drawer.draw(wrap=wrap)
