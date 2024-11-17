import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

# Kyte-Doolittle hydrophobicity scale
hydrophobicity_scale = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
    'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Function to parse FASTA input for protein sequence
def parse_fasta(fasta_text):
    fasta_io = StringIO(fasta_text)
    record = next(SeqIO.parse(fasta_io, "fasta"))
    return record

# Function to analyze protein sequence
def analyze_protein(sequence):
    protein_analysis = ProteinAnalysis(str(sequence))
    
    # Calculate key properties
    mol_weight = protein_analysis.molecular_weight()
    iso_point = protein_analysis.isoelectric_point()
    aa_count = protein_analysis.count_amino_acids()
    aromaticity = protein_analysis.aromaticity()
    instability_index = protein_analysis.instability_index()
    gravy = protein_analysis.gravy()
    helix, turn, sheet = protein_analysis.secondary_structure_fraction()
    
    # Hydropathy profile
    hydropathy_profile = [hydrophobicity_scale.get(aa, 0) for aa in sequence]
    
    return {
        "Molecular Weight (Da)": mol_weight,
        "Isoelectric Point": iso_point,
        "Aromaticity": aromaticity,
        "Instability Index": instability_index,
        "Hydrophobicity (GRAVY)": gravy,
        "Helix Content": helix,
        "Turn Content": turn,
        "Sheet Content": sheet,
        "Amino Acid Composition": aa_count,
        "Hydropathy Profile": hydropathy_profile
    }

# Function to plot hydropathy profile
def plot_hydropathy(hydropathy_profile):
    st.write("The hydropathy profile shows the hydrophobic and hydrophilic regions in the protein sequence. "
             "Peaks indicate hydrophobic regions, which may be membrane-associated, while valleys suggest hydrophilic regions, "
             "which are often exposed to the aqueous environment.")
    plt.figure(figsize=(10, 4))
    plt.plot(hydropathy_profile, color="blue")
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Hydropathy Index")
    plt.title("Hydropathy Profile")
    st.pyplot(plt)

# Function to plot secondary structure pie chart
def plot_secondary_structure_pie(helix, turn, sheet):
    st.write("The pie chart illustrates the predicted secondary structure composition of the protein. "
             "It shows the proportion of alpha-helices, beta-sheets, and turns, which are essential elements "
             "that contribute to the overall shape and function of the protein.")
    plt.figure(figsize=(6, 6))
    labels = ['Helix', 'Turn', 'Sheet']
    sizes = [helix, turn, sheet]
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.title("Secondary Structure Composition")
    st.pyplot(plt)

# Function to plot amino acid composition
def plot_aa_composition(aa_count):
    st.write("This bar chart displays the composition of each amino acid in the protein sequence. "
             "It helps identify which amino acids are most or least abundant, providing insights into "
             "the structural and functional properties of the protein.")
    plt.figure(figsize=(10, 5))
    plt.bar(aa_count.keys(), aa_count.values())
    plt.xlabel("Amino Acid")
    plt.ylabel("Count")
    plt.title("Amino Acid Composition")
    st.pyplot(plt)

# Educational Descriptions
descriptions = {
    "Molecular Weight (Da)": "Indicates the mass of the protein, essential for protein purification and analysis.",
    "Isoelectric Point": "The pH at which the protein has no net charge, affecting solubility and interaction with molecules.",
    "Aromaticity": "Represents the abundance of aromatic amino acids, often linked with protein stability.",
    "Instability Index": "Predicts stability in test conditions. Proteins with values above 40 are generally unstable.",
    "Hydrophobicity (GRAVY)": "The GRAVY (Grand Average of Hydropathy) value indicates hydrophobic or hydrophilic nature. Positive values suggest hydrophobicity, useful for membrane association predictions.",
    "Helix Content": "Shows the percentage of the protein likely to form alpha-helices, a common structural element.",
    "Turn Content": "Indicates regions that may form turns or bends, linking structural elements.",
    "Sheet Content": "Shows the proportion of beta-sheet structures, another key secondary structure element."
}

# App layout
st.title("Educational Protein Sequence Analyzer")

# Input FASTA sequence
fasta_input = st.text_area("Paste your protein sequence in FASTA format here", height=200)

# Submit button
if st.button("Analyze Sequence"):
    if fasta_input:
        try:
            # Parse the sequence
            record = parse_fasta(fasta_input)
            sequence = record.seq
            st.write(f"**Sequence ID**: {record.id}")
            st.write(f"**Sequence Description**: {record.description}")
            
            # Analyze the protein sequence
            analysis_results = analyze_protein(sequence)
            
            # Display results in a table with educational explanations
            summary_data = {
                "Property": [
                    "Molecular Weight (Da)", "Isoelectric Point", "Aromaticity",
                    "Instability Index", "Hydrophobicity (GRAVY)", "Helix Content",
                    "Turn Content", "Sheet Content"
                ],
                "Value": [
                    analysis_results["Molecular Weight (Da)"],
                    analysis_results["Isoelectric Point"],
                    analysis_results["Aromaticity"],
                    analysis_results["Instability Index"],
                    analysis_results["Hydrophobicity (GRAVY)"],
                    f"{analysis_results['Helix Content']:.2f}",
                    f"{analysis_results['Turn Content']:.2f}",
                    f"{analysis_results['Sheet Content']:.2f}"
                ],
                "Description": [
                    descriptions["Molecular Weight (Da)"],
                    descriptions["Isoelectric Point"],
                    descriptions["Aromaticity"],
                    descriptions["Instability Index"],
                    descriptions["Hydrophobicity (GRAVY)"],
                    descriptions["Helix Content"],
                    descriptions["Turn Content"],
                    descriptions["Sheet Content"]
                ]
            }
            summary_df = pd.DataFrame(summary_data)
            st.table(summary_df)
            
            # Display amino acid composition in a separate table
            aa_comp_df = pd.DataFrame(
                list(analysis_results["Amino Acid Composition"].items()),
                columns=["Amino Acid", "Count"]
            )
            st.subheader("Amino Acid Composition")
            st.table(aa_comp_df)
            
            # Plot amino acid composition
            st.subheader("Amino Acid Composition Plot")
            plot_aa_composition(analysis_results["Amino Acid Composition"])

            # Plot hydropathy profile
            st.subheader("Hydropathy Profile")
            plot_hydropathy(analysis_results["Hydropathy Profile"])

            # Plot secondary structure composition pie chart
            st.subheader("Secondary Structure Composition")
            plot_secondary_structure_pie(
                analysis_results["Helix Content"],
                analysis_results["Turn Content"],
                analysis_results["Sheet Content"]
            )

        except Exception as e:
            st.error(f"An error occurred: {e}")
    else:
        st.warning("Please enter a protein sequence in FASTA format.")
