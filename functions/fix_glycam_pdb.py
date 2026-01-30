from moleculekit.molecule import Molecule
import numpy as np


mol = Molecule('../glycosylation/3s5y_amber.pdb')

def safe_renumber(mol, chain, start_resid):
    """
    Renumbers a chain sequentially starting from start_resid 
    without causing ID collisions.
    """
    selection = f'chain {chain} and protein'
    original_resids = mol.get('resid', sel=selection)
    
    unique_old_resids, first_indices = np.unique(original_resids, return_index=True)
    unique_old_resids_ordered = original_resids[np.sort(first_indices)]
    
    mapping = {old: start_resid + i for i, old in enumerate(unique_old_resids_ordered)}

    new_resids = np.array([mapping[old] for old in original_resids])

    mol.set('resid', new_resids, sel=selection)

# Chain A: target range 32-421
safe_renumber(mol, 'A', 32)

# Chain B: target range 32-422
safe_renumber(mol, 'B', 32)

# Verification
check_a = np.unique(mol.get('resid', 'chain A and protein'))
check_b = np.unique(mol.get('resid', 'chain B and protein'))

print(f'Chain A range: {check_a.min()} to {check_a.max()}')
print(f'Chain B range: {check_b.min()} to {check_b.max()}')

mol.write('../glycosylation/3s5y_amber_fixed.pdb')