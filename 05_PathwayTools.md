# Microbiome Analysis with Pathway-Tools Integration

This enhanced microbiome analysis script now includes Pathway-Tools integration for comprehensive metabolic pathway analysis.

## Features Added

### 1. Metabolic Pathway Analysis
- **Pathway Abundance Calculation**: Maps taxonomic abundances to metabolic pathways
- **Differential Pathway Analysis**: Identifies pathways that differ between case and control groups
- **Pathway Volcano Plots**: Visualizes significant pathway differences
- **Pathway Heatmaps**: Shows pathway abundance patterns across groups

### 2. Metabolic Network Analysis
- **Pathway Correlation Networks**: Creates networks based on pathway correlations
- **Network Metrics**: Calculates network topology metrics
- **Network Visualization**: Generates interactive network plots

### 3. Pathway Enrichment Analysis
- **Subsystem Enrichment**: Analyzes enrichment of metabolic subsystems
- **Functional Annotation**: Maps pathways to biological functions
- **Enrichment Visualization**: Creates enrichment plots

### 4. Pathway-Tools Integration
- **Installation Detection**: Automatically detects Pathway-Tools installation
- **Script Generation**: Creates custom Pathway-Tools analysis scripts
- **API Integration**: Uses Pathway-Tools Python API when available

## Installation

### Prerequisites
1. **Pathway-Tools**: Install Pathway-Tools from [SRI International](https://bioinformatics.ai.sri.com/ptools/)
2. **Python Dependencies**: Install required Python packages

```bash
# Install Python dependencies
pip install -r requirements.txt

# For Pathway-Tools Python API (if available)
# pip install ptools
```

### Pathway-Tools Setup
1. Download and install Pathway-Tools
2. Add Pathway-Tools to your system PATH
3. Verify installation:
```bash
pathway-tools --version
```

## Usage

### Basic Usage
```bash
python 1_additional_microbiome_analysis.py
```

### Advanced Usage with Pathway-Tools
The script will automatically:
1. Check for Pathway-Tools installation
2. Perform metabolic pathway analysis
3. Generate Pathway-Tools scripts for further analysis
4. Create comprehensive metabolic reports

## Output Files

### New Metabolic Analysis Outputs
- `metabolic_analysis/metabolic_pathway_abundances.csv`: Pathway abundance data
- `metabolic_analysis/metabolic_pathway_analysis.csv`: Differential pathway analysis results
- `metabolic_analysis/pathway_volcano_plot.png`: Volcano plot of pathway differences
- `metabolic_analysis/pathway_heatmap.png`: Heatmap of pathway abundances
- `metabolic_analysis/metabolic_network.png`: Metabolic pathway network
- `metabolic_analysis/metabolic_network_metrics.csv`: Network topology metrics
- `metabolic_analysis/subsystem_enrichment.png`: Subsystem enrichment plot
- `metabolic_analysis/subsystem_enrichment_analysis.csv`: Enrichment analysis results
- `metabolic_analysis/pathway_tools_analysis.py`: Generated Pathway-Tools script

### Pathway-Tools Script
The generated `pathway_tools_analysis.py` script can be customized for:
- Detailed metabolic pathway analysis
- Reaction network reconstruction
- Compound analysis
- Flux balance analysis

## Metabolic Pathways Analyzed

The script analyzes the following metabolic pathways:

1. **Glycolysis**: Energy metabolism
2. **TCA Cycle**: Central carbon metabolism
3. **Amino Acid Metabolism**: Protein metabolism
4. **Fatty Acid Metabolism**: Lipid metabolism
5. **Vitamin Synthesis**: Vitamin production
6. **Antibiotic Resistance**: Resistance mechanisms
7. **Short Chain Fatty Acid Production**: SCFA synthesis
8. **Bile Acid Metabolism**: Bile acid processing
9. **Neurotransmitter Metabolism**: Neurotransmitter synthesis
10. **Inflammation Modulation**: Inflammatory response regulation

## Metabolic Subsystems

The enrichment analysis covers:

1. **Energy Metabolism**: Glycolysis, fermentation, respiration
2. **Amino Acid Metabolism**: Biosynthesis, degradation, transport
3. **Lipid Metabolism**: Fatty acid synthesis, bile acid metabolism
4. **Vitamin Metabolism**: Vitamin B and K synthesis
5. **Neurotransmitter Metabolism**: Serotonin and GABA metabolism

## Customization

### Adding New Pathways
Edit the `metabolic_pathways` dictionary in `create_metabolic_network_dataframe()`:

```python
metabolic_pathways = {
    'Your_Pathway': ['Taxon1', 'Taxon2', 'Taxon3'],
    # Add more pathways...
}
```

### Adding New Subsystems
Edit the `metabolic_subsystems` dictionary in `generate_pathway_enrichment_analysis()`:

```python
metabolic_subsystems = {
    'Your_Subsystem': {
        'taxa': ['Taxon1', 'Taxon2'],
        'functions': ['function1', 'function2']
    },
    # Add more subsystems...
}
```

### Pathway-Tools Script Customization
The generated `pathway_tools_analysis.py` script can be customized for:
- Specific metabolic analyses
- Custom pathway databases
- Advanced network reconstructions

## Troubleshooting

### Pathway-Tools Not Found
If Pathway-Tools is not detected:
1. Ensure Pathway-Tools is installed
2. Add to system PATH
3. Verify with `pathway-tools --version`

### Python API Issues
If the Pathway-Tools Python API is not available:
- The script will use command-line interface
- All analyses will still work with simplified metabolic analysis
- Custom scripts can be generated for manual Pathway-Tools analysis

### Memory Issues
For large datasets:
- Reduce the number of pathways analyzed
- Use sampling for network analysis
- Increase system memory allocation

## Citation

If you use this enhanced analysis in your research, please cite:
- Original microbiome analysis methods
- Pathway-Tools (if used): Karp et al. (2019) The Pathway Tools Pathway Prediction Algorithm. Stand Genomic Sci.
- Additional references for specific analyses used

## Support

For issues related to:
- **Pathway-Tools**: Contact SRI International
- **Python Analysis**: Check the script documentation
- **Custom Analysis**: Modify the script functions as needed 