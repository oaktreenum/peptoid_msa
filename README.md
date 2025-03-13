# Peptoid MSA Viewer  

Peptoid MSA Viewer is a Streamlit-based web application designed for visualizing and analyzing Multiple Sequence Alignments (MSA) of peptoids. The tool allows users to input MSA data, apply custom color mappings, and export results for further analysis.  

![MSA Example](images/MSA%20Example.png)  

## Features  

- **FASTA Upload & Parsing** – Supports input of MSA data in FASTA format.  
- **Custom Color Mapping** – Assign user-defined colors to specific peptoid residues.  
- **Dynamic Visualization** – Automatically adjusts for varying sequence lengths.  
- **Export Options** – Save alignment visualizations as PNG or SVG files.  
- **User-Friendly Interface** – Interactive design built with Streamlit.  

## Usage  

### Online Deployment  

This application is hosted on Hugging Face Spaces. You can access it directly via the following link:  

[Peptoid MSA Viewer on Hugging Face Spaces](https://huggingface.co/spaces/agoak/peptoid-msa)  

### Local Installation  

To run the application locally, follow these steps:  

1. Clone the repository:  

   ```
   git clone <link_to_repo>  
   cd peptoid-msa
   ```

1. Install dependencies:  

   ```
   pip install -r requirements.txt
   ```

2. Start the Streamlit app:  

   ```
   streamlit run app.py
   ```

## File Formats  

- **Input:** FASTA formatted 3-letter code peptoid sequences pasted into the input box.  
- **Output:** PNG, SVG, or FASTA download.  

## Dependencies  

The application requires the following Python packages:  

- numpy 
- plotly
- streamlit 
- kaleido

All dependencies are listed in `requirements.txt` and will be installed automatically.  

## Author  

**Allon Goldberg**  
Research Assistant  
Flatiron Institute, Center for Computational Biology (CCB)  
Biomolecular Design Group  

March 13, 2025  

## Acknowledgments  

If you use this tool in your research, please cite appropriately.  
