# packages/chirality_detector.py

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

class ChiralityDetector:
    def analyze(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Molécula inválida"}
        chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        return {"smiles": smiles, "chiral_centers": chiral_centers}
