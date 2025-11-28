import sys
import os
import networkx as nx
import json

# Add backend directory to path so we can import app
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from backend.app.parser import GraphParser

def test_smiles():
    print("Testing SMILES...")
    smiles = "CCO" # Ethanol
    g = GraphParser.parse(smiles)
    assert isinstance(g, nx.Graph)
    assert len(g.nodes) == 3
    assert len(g.edges) == 2
    # Check labels
    labels = nx.get_node_attributes(g, 'label')
    print(f"SMILES Nodes: {labels}")
    assert labels[0] == 'C'
    assert labels[2] == 'O'
    print("SMILES Test Passed!")

def test_json():
    print("\nTesting JSON...")
    g_orig = nx.Graph()
    g_orig.add_edge(0, 1)
    data = nx.node_link_data(g_orig)
    json_str = json.dumps(data)
    g = GraphParser.parse(json_str)
    assert isinstance(g, nx.Graph)
    assert len(g.nodes) == 2
    assert len(g.edges) == 1
    print("JSON Test Passed!")

def test_matrix():
    print("\nTesting Adjacency Matrix...")
    matrix_str = "[[0, 1], [1, 0]]"
    g = GraphParser.parse(matrix_str)
    assert isinstance(g, nx.Graph)
    assert len(g.nodes) == 2
    assert len(g.edges) == 1
    print("Matrix Test Passed!")

def test_edgelist():
    print("\nTesting Edge List...")
    edgelist_str = "1 2\n2 3"
    g = GraphParser.parse(edgelist_str)
    assert isinstance(g, nx.Graph)
    assert len(g.nodes) == 3
    assert len(g.edges) == 2
    print("Edge List Test Passed!")

if __name__ == "__main__":
    try:
        test_smiles()
        test_json()
        test_matrix()
        test_edgelist()
        print("\nALL TESTS PASSED")
    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        sys.exit(1)
