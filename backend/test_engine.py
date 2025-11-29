import sys
import os
import networkx as nx

# Add backend directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from backend.app.engine import analyze_graph

def test_planar_graph():
    print("Testing Planar Graph (K4)...")
    G = nx.complete_graph(4)
    result = analyze_graph(G)
    
    assert result['is_planar'] is True
    assert len(result['nodes']) == 4
    assert len(result['edges']) == 6
    
    # Check that no edges are marked as conflict
    for edge in result['edges']:
        assert edge['is_conflict'] is False
        
    print("Planar Graph Test Passed!")

def test_non_planar_graph():
    print("\nTesting Non-Planar Graph (K5)...")
    G = nx.complete_graph(5)
    result = analyze_graph(G)
    
    assert result['is_planar'] is False
    assert len(result['nodes']) == 5
    assert len(result['edges']) == 10
    
    # Check that some edges are marked as conflict
    conflict_count = sum(1 for edge in result['edges'] if edge['is_conflict'])
    print(f"Conflict Edges Found: {conflict_count}")
    assert conflict_count > 0
    
    print("Non-Planar Graph Test Passed!")

if __name__ == "__main__":
    try:
        test_planar_graph()
        test_non_planar_graph()
        print("\nALL TESTS PASSED")
    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        sys.exit(1)
