# RNA-Seq Pseudoalignment using De Bruijn Graphs

## Overview

This project implements a pseudoalignment algorithm for RNA-Seq data using colored De Bruijn graphs. The approach allows for efficient mapping of sequencing reads to transcript isoforms without performing full alignment, significantly reducing computational requirements.

## Key Features

- **Constructs a colored De Bruijn graph** from transcript isoforms  
- **Handles reverse complement reads** automatically  
- **Processes reads with ambiguous bases** ('N')  
- **Generates equivalence classes** for quantification  
- **Outputs results** in a structured CSV format  
