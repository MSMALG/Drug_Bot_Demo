import pandas as pd
import streamlit as st

file_path = 'drug_nomenclature_and_modifications.xlsx'  
df = pd.read_excel(file_path)

#Function to fetch molecule image from PubChem
def fetch_pubchem_image(smiles, width=400, height=400):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/png"
    #Create the request URL
    url = f"{base_url}?smiles={smiles}&width={width}&height={height}"
    return url

st.title("Drug Recommendation and Side Effect Solution Bot")

target_receptor = st.text_input("Enter the target receptor (e.g., β1-selective adrenergic blocker):")

if target_receptor:
    filtered_df = df[df['Target Receptor'].str.contains(target_receptor, case=False, na=False)]
    
    if not filtered_df.empty:
        drug_name = filtered_df.iloc[0]['Drug Nomenclature']
        smiles = filtered_df.iloc[0]['Chemical Nomenclature (Molecular Input Line Entry System)']
        
        st.subheader(f"Recommended Drug: {drug_name}")
        
        #Display the molecular structure 
        if smiles:
            pubchem_url = fetch_pubchem_image(smiles)
            st.image(pubchem_url, caption=f"Molecular Structure of {drug_name}")
        else:
            st.error("SMILES string not available for this molecule.")

        if st.button("Check for Side Effects"):
            side_effects = filtered_df['Side Effect'].tolist()
            modifications = filtered_df['Modification'].tolist()
            modified_smiles = filtered_df['Molecular Input Line Entry System (Modified)'].tolist()
            
            for side_effect, modification, mod_smiles in zip(side_effects, modifications, modified_smiles):
                st.subheader(f"Side Effect: {side_effect}")
                st.text(f"Solution: {modification}")
                
                #Display the modified molecular structure 
                if mod_smiles:
                    pubchem_mod_url = fetch_pubchem_image(mod_smiles)
                    st.image(pubchem_mod_url, caption="Modified Molecular Structure")
                else:
                    st.error(f"Modified SMILES not available for the side effect: {side_effect}.")
    else:
        st.error("No drugs found for the specified target receptor.")
