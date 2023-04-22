import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


# Set the app title
st.title("Molecule Similarity Calculator")
    
# Get the SMILES notations from the user
smiles1 = st.text_input("Enter the SMILES notation for molecule 1:")
smiles2 = st.text_input("Enter the SMILES notation for molecule 2:")
    
# Generate RDKit molecules from the SMILES strings
mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)
    
# Calculate the ECFP fingerprints for each molecule
fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    
# Calculate the Tanimoto similarity index between the two fingerprints
similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
# Display the similarity index to the user
st.write(f"The Tanimoto similarity index between the two molecules is: {similarity:.2f}")
