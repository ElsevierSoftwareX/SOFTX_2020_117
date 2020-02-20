"""Model node graph visualisation with Graphviz.

Based on node_graph project by Sebastian Steinlechner: https://github.com/sestei/node_graph/
"""

import tempfile
from ...components import laser, space

class NodeGraph(object):
    """Pykat node graph plotter."""
    def __init__(self, kat):
        self.kat = kat

    def view_pdf(self):
        """View the graph as a PDF"""
        return self.node_graph().view(directory=tempfile.gettempdir(), cleanup=True)

    def _repr_svg_(self):
        """Graphviz rendering for Jupyter notebooks."""
        return self.node_graph()._repr_svg_()

    def node_graph(self, engine="fdp", node_style="filled", node_font="Helvetica",
                   node_font_size=10, node_width=0.2, node_shape="point", node_color="red",
                   node_dump_color="black", edge_arrowhead="dot", graph_font="Helvetica",
                   graph_font_size=8, component_color="LightSkyBlue", laser_color="Orange",
                   space_color="MediumSlateBlue", detector_color="YellowGreen"):
        """Create Graphviz node graph"""
        from graphviz import Digraph
        graph = Digraph(engine=engine)
        graph.attr("node", style=node_style, fontname=node_font, fontsize=str(node_font_size))
        graph.attr("edge", arrowhead=edge_arrowhead)
        graph.attr("graph", fontname=graph_font,
                   fontsize=str(graph_font_size))

        def add_kat_node(from_id, node):
            node_id = str(node.id)
            color = node_dump_color if node.name == 'dump' else node_color
            graph.node(node_id, label='', xlabel=node.name, width=str(node_width),
                       shape=node_shape, color=color)
            graph.edge(from_id, node_id)

        for comp in self.kat.components.values():
            if isinstance(comp, laser):
                attr = {'fillcolor': laser_color, 'shape': 'box'}
            elif isinstance(comp, space):
                attr = {'fillcolor': space_color, 'shape': 'diamond'}
            else:
                attr = {'fillcolor': component_color, 'shape': 'box'}
            graph.node(comp.name, **attr)
            for node in comp.nodes:
                add_kat_node(comp.name, node)

        for det in self.kat.detectors.values():
            if len(det._nodes) > 0:
                graph.node(det.name, fillcolor=detector_color, shape='ellipse')
                for node in det._nodes:
                    add_kat_node(det.name, node)

        return graph