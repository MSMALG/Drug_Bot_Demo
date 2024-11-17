import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

file_path = 'final_drugs.xlsx'
df = pd.read_excel(file_path)

st.title("Drug Recommendation and Side Effect Solution Bot")

if "show_user_guide" not in st.session_state:
    st.session_state["show_user_guide"] = False

if st.button("User Guide"):
    st.session_state["show_user_guide"] = not st.session_state["show_user_guide"]

if st.session_state["show_user_guide"]:
    st.markdown("""
     ## Instructions on How to Use the App:
    Welcome to our innovative prototype for drug discovery and design! 
    This app allows you to explore drugs for specific receptors and view modifications to address their side effects.

    Here's a quick guide to help you get started:

    1. **Search for Molecules**: Enter a receptor of interest in the search bar. 
       Currently, you can search for drugs targeting the following receptors:
       - β1-selective adrenergic blocker
       - β1-adrenergic agonist
       - β2-adrenergic antagonist
       - β2-adrenergic agonist
       - α1-adrenergic antagonist

    2. **Explore Side Effects**: Once a drug is displayed, you can type a query such as "What are the side effects?" 
       to view a list of side effects associated with the drug.

    3. **View Modifications**: Select a specific side effect to see the proposed solution 
       and the molecular structure of the modified drug.
    """)
    st.divider()  

# User inputs the target receptor
target_receptor = st.text_input("Enter the target receptor (e.g., β1-selective adrenergic blocker):")

if target_receptor:
    # Query the dataset for the given receptor
    filtered_df = df[df['Target Receptor'].str.contains(target_receptor, case=False, na=False)]
    
    if not filtered_df.empty:
        drug_name = filtered_df.iloc[0]['Drug Nomenclature']
        smiles = filtered_df.iloc[0]['Chemical Nomenclature (Molecular Input Line Entry System)']
        
        # Display the recommended drug
        st.subheader(f"Recommended Drug: {drug_name}")
        
        # Validate and display the molecular structure
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            st.image(img, caption=f"Molecular Structure of {drug_name}")
        else:
            st.error("Invalid SMILES string. Unable to generate molecular structure.")

        # Query box for user input
        user_query = st.text_input("Any other query?")

        if user_query:
            # Check if the query contains the words "side effects"
            if "side effects" in user_query.lower():
                # Extract side effects for the specific drug
                side_effects = filtered_df['Side Effect'].tolist()
                modifications = filtered_df['Modification'].tolist()
                modified_smiles = filtered_df['Molecular Input Line Entry System (Modified)'].tolist()
                
                if side_effects:
                    st.subheader(f"Side Effects for {drug_name}")
                    
                    # Display each side effect as an expandable section
                    for i, (side_effect, modification, mod_smiles) in enumerate(zip(side_effects, modifications, modified_smiles)):
                        with st.expander(f"Side Effect {i + 1}: {side_effect}"):
                            st.text(f"Solution: {modification}")
                            
                            # Validate and visualize modified molecule
                            mod_mol = Chem.MolFromSmiles(mod_smiles)
                            if mod_mol:
                                mod_img = Draw.MolToImage(mod_mol)
                                st.image(mod_img, caption="Modified Molecular Structure")
                            else:
                                st.error(f"Invalid SMILES string for modified molecule of {side_effect}.")
                else:
                    st.error(f"No side effects listed for {drug_name}.")
            else:
                st.info("Please specify a query about side effects to see relevant results.")
    else:
        st.error("No drugs found for the specified target receptor.")
