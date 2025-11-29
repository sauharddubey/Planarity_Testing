import networkx as nx
import json
import numpy as np
from rdkit import Chem
import re

class GraphParser:
    @staticmethod
    def parse(input_str: str) -> nx.Graph:
        input_str = input_str.strip()
        
        # 1. JSON Detection
        if input_str.startswith("{"):
            try:
                data = json.loads(input_str)
                # Check for "nodes" and "edges" format
                if "nodes" in data and "edges" in data:
                    g = nx.Graph()
                    for node in data["nodes"]:
                        # Add node with attributes (x, y, etc.)
                        g.add_node(node["id"], **node)
                    for edge in data["edges"]:
                        g.add_edge(edge["source"], edge["target"])
                    return g
                
                # Fallback to standard node_link_graph
                return nx.node_link_graph(data)
            except (json.JSONDecodeError, nx.NetworkXError, KeyError):
                pass # Fallthrough

        # 2. Adjacency Matrix Detection
        if input_str.startswith("[["):
            try:
                matrix = json.loads(input_str)
                np_matrix = np.array(matrix)
                return nx.from_numpy_array(np_matrix)
            except (json.JSONDecodeError, ValueError):
                pass # Fallthrough

        # 3. SMILES Detection
        # Heuristic: Contains typical chemical symbols, no spaces, not an edge list
        # Common SMILES chars: C, N, O, F, Cl, Br, I, S, P, =, #, (, ), [, ], +, -, @, /
        # Edge list usually has spaces. SMILES usually doesn't (except maybe isomeric SMILES which we can handle if needed, but standard ones don't)
        if " " not in input_str and re.search(r"[CNOcno]", input_str):
            try:
                mol = Chem.MolFromSmiles(input_str)
                if mol:
                    g = nx.Graph()
                    for atom in mol.GetAtoms():
                        g.add_node(atom.GetIdx(), label=atom.GetSymbol())
                    for bond in mol.GetBonds():
                        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    return g
            except Exception:
                pass # Fallthrough

        # 4. Edge List Fallback
        try:
            # nx.parse_edgelist expects lines.
            # If input is "1 2", it's one line.
            # If input is "1 2\n2 3", it's multiple.
            return nx.parse_edgelist(input_str.splitlines())
        except Exception as e:
            raise ValueError(f"Could not parse input as any known format. Error: {e}")
