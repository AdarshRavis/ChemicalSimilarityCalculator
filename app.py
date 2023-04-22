import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Set the app title
st.title("Molecular Similarity Calculator")

st.write("Adarsh Ravishankar, MD")
st.write("Departments of Medicine and Dermatology, University of Minnesota")

# Get the SMILES notations from the user
smiles1 = st.text_input("Enter the SMILES notation for molecule 1:")
smiles2 = st.text_input("Enter the SMILES notation for molecule 2:")

# Define the "Calculate" button
if st.button("Calculate"):
    # Generate RDKit molecules from the SMILES strings
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
    except:
        st.error("Invalid SMILES notation entered. Please enter valid SMILES notation and try again.")
        st.stop()
        
    # Calculate the ECFP fingerprints for each molecule
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    
    # Calculate the Tanimoto similarity index between the two fingerprints
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    # Display the similarity index to the user
    st.markdown(f"**The Tanimoto similarity index between the two molecules is: {similarity:.2f}**")
