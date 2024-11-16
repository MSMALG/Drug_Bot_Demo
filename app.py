import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

file_path = 'drug_nomenclature_and_modifications.xlsx'
df = pd.read_excel(file_path)

st.title("Drug Recommendation and Side Effect Solution Bot")

target_receptor = st.text_input("Enter the target receptor (e.g., β1-selective adrenergic blocker):")

if target_receptor:
    #Query the dataset for the given receptor
    filtered_df = df[df['Target Receptor'].str.contains(target_receptor, case=False, na=False)]
    
    if not filtered_df.empty:
        drug_name = filtered_df.iloc[0]['Drug Nomenclature']
        smiles = filtered_df.iloc[0]['Chemical Nomenclature (Molecular Input Line Entry System)']
        
        st.subheader(f"Recommended Drug: {drug_name}")
        
        #Validate and display the molecular structure
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            st.image(img, caption=f"Molecular Structure of {drug_name}")
        else:
            st.error("Invalid SMILES string. Unable to generate molecular structure.")

        if st.button("Check for Side Effects"):
            side_effects = filtered_df['Side Effect'].tolist()
            modifications = filtered_df['Modification'].tolist()
            modified_smiles = filtered_df['Molecular Input Line Entry System (Modified)'].tolist()
            
            for side_effect, modification, mod_smiles in zip(side_effects, modifications, modified_smiles):
                st.subheader(f"Side Effect: {side_effect}")
                st.text(f"Solution: {modification}")
                
                #Validate and visualize Modified Molecule
                mod_mol = Chem.MolFromSmiles(mod_smiles)
                if mod_mol:
                    mod_img = Draw.MolToImage(mod_mol)
                    st.image(mod_img, caption="Modified Molecular Structure")
                else:
                    st.error(f"Invalid SMILES string for modified molecule of {side_effect}.")
    else:
        st.error("No drugs found for the specified target receptor.")
