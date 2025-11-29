import logging
import json
from typing import Dict, Any, Optional, Union, List
from concurrent.futures import Executor
import asyncio
from app.worker import process_graph_task

# RDKit imports
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

import os
from dotenv import load_dotenv
import google.generativeai as genai

# Load environment variables
load_dotenv()

# Configure Gemini
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)
else:
    logger.warning("GEMINI_API_KEY not found in environment variables.")

logger = logging.getLogger(__name__)



def smiles_to_edge_list(mol) -> str:
    """Converts an RDKit molecule to a string edge list (0-indexed)."""
    edges = []
    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        edges.append(f"{u} {v}")
    return "\n".join(edges)

def analyze_molecule(smiles: str) -> Dict[str, Any]:
    """
    Analyzes a molecule given its SMILES string using RDKit.
    Also runs Planarity Testing on the molecular graph.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"valid": False, "error": "Invalid SMILES string"}

        # Calculate properties
        mol_wt = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Lipinski's Rule of 5 check
        lipinski_violations = 0
        if mol_wt > 500: lipinski_violations += 1
        if logp > 5: lipinski_violations += 1
        if hbd > 5: lipinski_violations += 1
        if hba > 10: lipinski_violations += 1

        # Generate 3D coordinates (keep for potential future use or fallback)
        mol_3d = Chem.AddHs(mol)
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        
        # Extract atoms and bonds for 3D rendering
        atoms_3d = []
        conf_3d = mol_3d.GetConformer()
        for atom in mol_3d.GetAtoms():
            pos = conf_3d.GetAtomPosition(atom.GetIdx())
            atoms_3d.append({
                "symbol": atom.GetSymbol(),
                "x": pos.x,
                "y": pos.y,
                "z": pos.z,
                "color": get_atom_color(atom.GetSymbol())
            })
            
        bonds_3d = []
        for bond in mol_3d.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            start_pos = conf_3d.GetAtomPosition(start_idx)
            end_pos = conf_3d.GetAtomPosition(end_idx)
            bonds_3d.append({
                "start": {"x": start_pos.x, "y": start_pos.y, "z": start_pos.z},
                "end": {"x": end_pos.x, "y": end_pos.y, "z": end_pos.z},
                "order": bond.GetBondTypeAsDouble()
            })

        # --- Generate 2D Coordinates for D3 ---
        mol_2d = Chem.Mol(mol) # Copy
        AllChem.Compute2DCoords(mol_2d)
        
        atoms_2d = []
        conf_2d = mol_2d.GetConformer()
        for atom in mol_2d.GetAtoms():
            pos = conf_2d.GetAtomPosition(atom.GetIdx())
            atoms_2d.append({
                "id": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "x": pos.x, # RDKit 2D coords are usually small floats, frontend will scale
                "y": pos.y,
                "color": get_atom_color(atom.GetSymbol())
            })

        bonds_2d = []
        for bond in mol_2d.GetBonds():
            bonds_2d.append({
                "source": bond.GetBeginAtomIdx(),
                "target": bond.GetEndAtomIdx(),
                "order": bond.GetBondTypeAsDouble()
            })

        # --- Planarity Check ---

        # --- Planarity Check ---
        edge_list_str = smiles_to_edge_list(mol)
        # We run this synchronously here because we are already in a thread (run_in_executor)
        # But process_graph_task is a regular function, so it's fine.
        planarity_response = process_graph_task(edge_list_str)
        
        planarity_result = {}
        if planarity_response.get("status") == "success":
            planarity_result = planarity_response.get("data", {})
        else:
            logger.warning(f"Planarity check failed: {planarity_response.get('message')}")
            planarity_result = {"is_planar": True, "error": "Check failed"} # Default to safe assumption or handle error

        return {
            "valid": True,
            "smiles": smiles,
            "properties": {
                "molecular_weight": round(mol_wt, 2),
                "logp": round(logp, 2),
                "hbd": hbd,
                "hba": hba,
                "tpsa": round(tpsa, 2),
                "rotatable_bonds": rotatable_bonds,
            },
            "lipinski": {
                "violations": lipinski_violations,
                "passes": lipinski_violations <= 1
            },
            "structure_3d": {
                "atoms": atoms_3d,
                "bonds": bonds_3d
            },
            "structure_2d": {
                "nodes": atoms_2d,
                "links": bonds_2d
            },
            "planarity": planarity_result
        }
    except Exception as e:
        logger.error(f"Error analyzing molecule: {e}")
        return {"valid": False, "error": str(e)}

def get_atom_color(symbol: str) -> str:
    # CPK coloring convention
    colors = {
        "C": "#909090", # Carbon - Grey
        "H": "#FFFFFF", # Hydrogen - White
        "O": "#FF0D0D", # Oxygen - Red
        "N": "#3050F8", # Nitrogen - Blue
        "S": "#FFFF30", # Sulfur - Yellow
        "P": "#FFA500", # Phosphorus - Orange
        "F": "#90E050", # Fluorine - Green
        "Cl": "#1FF01F", # Chlorine - Green
        "Br": "#A62929", # Bromine - Dark Red
        "I": "#940094", # Iodine - Violet
    }
    return colors.get(symbol, "#FF00FF") # Default Magenta

async def generate_llm_response(analysis: Union[Optional[Dict[str, Any]], List[Dict[str, Any]]], user_query: str) -> str:
    """
    Generates a response using Google Gemini based on the analysis and user query.
    """
    if not GEMINI_API_KEY:
        return "Error: GEMINI_API_KEY not configured."

    try:
        model = genai.GenerativeModel('gemini-2.0-flash')
        
        system_prompt = (
            "You are an expert Pharmaceutical Scientist and Chemist assistant designed to aid in drug discovery. "
            "Your role is to analyze molecular properties, explain structure-activity relationships (SAR), and suggest optimizations.\n"
            "Guidelines:\n"
            "- Be scientific, precise, and professional.\n"
            "- If a molecule is analyzed, interpret its properties (MW, LogP, Lipinski rules) in the context of drug-likeness.\n"
            "- If multiple molecules are provided, COMPARE them. Highlight the best candidates based on drug-likeness and the user's query.\n"
            "- Highlight potential risks (e.g., toxicity, poor solubility) or benefits.\n"
            "- **Topological Analysis**: The system also checks if the molecular graph is PLANAR. Mention this. If non-planar (contains K5 or K3,3 subgraphs), explain that this implies a complex, potentially caged or cross-linked structure.\n"
            "- If no molecule is provided, answer general chemistry questions or ask for a SMILES string.\n"
            "- Keep responses concise and actionable."
        )

        prompt_parts = [system_prompt]

        if analysis:
            if isinstance(analysis, list):
                # Batch Analysis Context
                context = "Batch Analysis of Molecules:\n"
                for i, res in enumerate(analysis):
                    if res.get("valid"):
                        props = res['properties']
                        lip = res['lipinski']
                        context += (
                            f"\nMolecule {i+1} (SMILES: {res['smiles']}):\n"
                            f"- MW: {props['molecular_weight']}, LogP: {props['logp']}, "
                            f"HBD: {props['hbd']}, HBA: {props['hba']}, TPSA: {props['tpsa']}\n"
                            f"- Lipinski: {'Pass' if lip['passes'] else 'Fail'} ({lip['violations']} violations)\n"
                            f"- Planarity: {'Planar' if res['planarity']['is_planar'] else 'Non-Planar'}\n"
                        )
                        if not res['planarity']['is_planar']:
                            context += f"  - Kuratowski Subgraphs: K5={res['planarity'].get('k5_count', 0)}, K3,3={res['planarity'].get('k33_count', 0)}\n"
                    else:
                        context += f"\nMolecule {i+1} (SMILES: {res.get('smiles', 'Unknown')}): Invalid or Error ({res.get('error')})\n"
                prompt_parts.append(f"Context:\n{context}\n\nUser Query: {user_query}")
            elif isinstance(analysis, dict) and analysis.get("valid"):
                # Single Molecule Context
                props = analysis['properties']
                lip = analysis['lipinski']
                context = (
                    f"Current Molecule Analysis (SMILES: {analysis['smiles']}):\n"
                    f"- Molecular Weight: {props['molecular_weight']} g/mol\n"
                    f"- LogP: {props['logp']}\n"
                    f"- H-Bond Donors: {props['hbd']}\n"
                    f"- H-Bond Acceptors: {props['hba']}\n"
                    f"- TPSA: {props['tpsa']}\n"
                    f"- Lipinski Rule of 5: {'Passes' if lip['passes'] else 'Fails'} ({lip['violations']} violations)\n"
                    f"- Graph Topology: {'Planar' if analysis['planarity']['is_planar'] else 'Non-Planar'}\n"
                )
                if not analysis['planarity']['is_planar']:
                    context += f"- Kuratowski Subgraphs: K5={analysis['planarity'].get('k5_count', 0)}, K3,3={analysis['planarity'].get('k33_count', 0)}\n"
                prompt_parts.append(f"Context:\n{context}\n\nUser Query: {user_query}")
        else:
            prompt_parts.append(f"User Query: {user_query}")

        response = await model.generate_content_async(prompt_parts)
        
        return response.text
    except Exception as e:
        logger.error(f"Error calling Gemini: {e}")
        return f"Error generating response: {str(e)}"
