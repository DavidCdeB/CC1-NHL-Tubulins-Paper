from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename="complex_start_outer.pdb")

# Do NOT add missing residues (keeps your numbering stable)
fixer.findMissingResidues()
fixer.missingResidues = {}

# Add missing heavy atoms within existing residues
fixer.findMissingAtoms()
fixer.addMissingAtoms()

with open("complex_fixed_atoms.pdb", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
