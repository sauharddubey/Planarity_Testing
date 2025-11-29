import networkx as nx

def analyze_graph(G: nx.Graph) -> dict:
    """
    Analyzes the graph for planarity and generates a layout.
    Returns a dictionary with nodes, edges, and planarity status.
    """
    is_planar, certificate = nx.check_planarity(G, counterexample=True)
    
    conflict_edges = set()
    if not is_planar:
        # certificate is the counterexample subgraph (Kuratowski subgraph)
        for u, v in certificate.edges():
            conflict_edges.add(tuple(sorted((u, v))))
    
    # Determine layout
    # 1. Check if graph already has coordinates (from input JSON)
    has_coords = all('x' in G.nodes[n] and 'y' in G.nodes[n] for n in G.nodes())
    
    if has_coords:
        # Use existing coordinates
        pos = {n: (G.nodes[n]['x'], G.nodes[n]['y']) for n in G.nodes()}
    else:
        # 2. Generate layout
        if not is_planar:
            # Use spring layout for non-planar graphs
            pos = nx.spring_layout(G, seed=42)
        else:
            # Use planar layout for planar graphs
            try:
                pos = nx.planar_layout(G)
            except nx.NetworkXException:
                pos = nx.spring_layout(G, seed=42)

    nodes = []
    for node in G.nodes():
        # pos[node] is a numpy array or list [x, y]
        x, y = pos[node]
        # Get label if exists, else use node ID
        label = str(G.nodes[node].get('label', node))
        nodes.append({
            "id": node,
            "x": float(x),
            "y": float(y),
            "label": label
        })
        
    edges = []
    for u, v in G.edges():
        # Check if this edge is in the conflict set
        is_conflict = tuple(sorted((u, v))) in conflict_edges
        edges.append({
            "source": u,
            "target": v,
            "is_conflict": is_conflict
        })
        
    return {
        "is_planar": is_planar,
        "nodes": nodes,
        "edges": edges
    }
