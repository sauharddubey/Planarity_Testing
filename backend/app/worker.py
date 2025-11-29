from app.parser import GraphParser
from app.engine import analyze_graph
import networkx as nx

def process_graph_task(input_str: str):
    try:
        # Parse the graph
        graph = GraphParser.parse(input_str)
        
        # Analyze the graph
        result = analyze_graph(graph)
        
        return {"status": "success", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}
