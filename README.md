The BSH project aims to identify and map protein sequences from sequencing data. Additionally, it focuses on specific enzymes listed in the KEGG and Rhea databases, including proteins that have not yet been annotated.

## Analysis Workflow

1. **Prepare datasets:**
   - Representative protein sequences from clusters
   - An indexed protein database
   - Enzyme sequences

2. **Run BLAST** to search the protein database using the representative protein sequences.

3. **Identify enzymes** in the BLAST results.

4. **Aggregate enzymes with the same function**, sum their read counts, and determine whether they are enriched in the dataset.
