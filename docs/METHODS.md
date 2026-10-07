# Methods

The input parser preserves IDs, source line numbers and observed calls. Calls are sorted to normalize allele order. Coordinates receive a consensus only when observed calls agree; missing calls do not overwrite observations. Mixed indel/base encodings and conflicting calls remain in separate audits.

Seven archived Ensembl / NCBI coordinate anchors support GRCh37 for the original export. NCBI SPDI coordinates are converted from zero-based to one-based. This does not identify the chip, provider or strand convention. The workflow does not lift coordinates between assemblies.

The full public 1000 Genomes Phase 3 PGEN/PVAR/PSAM files are pinned by length and SHA-256. Actual PSAM IDs and both label levels are checked against the official panel. Every selected variant retains its PVAR index, reference/alternate alleles, sample source ID and source line. At each coordinate, exactly one eligible biallelic SNP must be available. Palindromic A/T and C/G pairs are excluded; a unique reverse complement is accepted and recorded. Ambiguous or incompatible sites are excluded.

Global references comprise all 2,504 individuals; the East Asian fit uses only 504 EAS individuals and labels CHB, CHS, JPT, CDX and KHV. Each scope selects complete-call variants independently. MAF thresholds are 0.01 globally and 0.05 within East Asia. A greedy Pearson dosage-correlation filter follows PVAR order, excluding a site when r² exceeds 0.2 against a retained site on the same chromosome within 500 kb. This is an explicitly implemented method, not a claim of identical PLINK pruning defaults. Population structure can contribute to correlation.

Reference genotypes are standardized as `(g − 2p) / sqrt(2p(1 − p))`. The target never changes frequencies, selection, LD or the reference PCA. Four components are fitted, and the target is transformed through the fixed model. The global fit uses randomized SVD with seed 20261006, 20 oversamples and seven power iterations. The East Asian fit uses full SVD to reduce approximation effects on small components. Float32 computation is independently reconstructed and checked in float64.

Distances use PC1–PC4 scaled by reference PC standard deviations. A two-PC comparison is also recorded. Neither distance is converted to an ancestry fraction or confidence probability.

The 80/20 diagnostic uses stratified reference splitting. Feature selection, frequencies, LD and PCA are refitted on training references only; held-out individuals are compared with training centroids. Per-group counts and agreement are retained. Leave-one-chromosome-out refits remove each autosome from the main selected marker set without repeating selection. Its counts describe model sensitivity, not a confidence interval.

Low holdout agreement for an adjacent reference group and changes in the target's nearest center limit fine-scale interpretation even when broad global placement is stable. Reference ascertainment, unknown chip/strand metadata and limited marker coverage remain relevant. This workflow performs no admixture, ancient ancestry, haplogroup or clinical inference.
