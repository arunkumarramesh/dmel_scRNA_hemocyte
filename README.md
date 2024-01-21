# Drosophila larval hemocyte single cell RNA-seq
Scripts for "Constitutive activation of cellular immunity underlies the evolution of resistance to infection in Drosophila"

Alexandre B Leitão*, Ramesh Arunkumar*, Jonathan P Day,Emma M Geldman,Ismaël Morin-Poulard, Michèle Crozatier, Francis M Jiggins 
https://doi.org/10.7554/eLife.59095

# Description
Goal was to identify cell types in Drosophila larval hemocytes (blood cells) from 10X single cell RNA-sequencing using Seurat

# List of scripts

singlecell_dmel_hemocytes.R - Single cell analysis

DoHeatmap2.R - Modified Seurat v3 DoHeatmap function for plotting

compare_datasets.R - Comparing three Drosophilia single cell hemocyte datasets

app.R - Shiny app for single cell gene expression browser

# Overview of Seurat clustering pipeline 
also in https://doi.org/10.7554/eLife.59095

![lax_59095_elife-59095-fig3-figsupp6-v2 tif](https://github.com/arunkumarramesh/dmel_scRNA_hemocyte/assets/23363383/7558021d-a96e-40fc-bb7e-738795f489c8)

# Gene expression browser
https://arunkuma.shinyapps.io/singlecell_app

Goal is to identify genes that are constitutively activated following long-term parasite exposure and those that are transiently induced by parasite infection

- Input gene name as Flybase gene ID, which can be obtained from https://flybase.org/
- PLASM1-2: plasmatocytes, LAM1-2: immature lamellocytes, LAM3: mature lamellocytes, CC: crystal cells, MET: metabolism, AMP: anti-microbial peptide
- Selection: 26 generations of selection under high parasite (Leptopilina boulardi) pressure
- No Selection: 26 generations under standard conditions without parasite exposure
- Infection: Infection with parasitic wasp
- No Infection: No infection control
- Expression level: Log-normalised read counts from Seurat V3
  
Full data available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148826

Please post questions at https://github.com/arunkumarramesh/dmel_scRNA_hemocyte


<img width="1391" alt="Screenshot 2024-01-20 at 21 53 57" src="https://github.com/arunkumarramesh/dmel_scRNA_hemocyte/assets/23363383/e7e1b78e-22db-42a7-b92a-c70bb0cf9505">




