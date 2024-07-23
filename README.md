# Phasepaint

`phasepaint` is a tool for correcting phasing errors after ancestry painting of genomic data.

Ancestry painting method methods (such as [Loter]()) usually require phased haplotypes. For diploid, phase is usually inferred using tools like [whatshap]() and [shapeit](), but the inference is imperfect. "Phase switch errors" occur when the phasing tool incorrectly swaps the maternal and paternal haplotypes. When combined with ancestry painting, phase switch errors can appear as breaks in ancestry blocks.

`phasepaint` aims to correct phase switch errors after ancestry painting. It uses a simple heuristic algorithm based on the assumption that the optimal phasing is that which maximises similarity between haplotypes in the dataset.

No installation is required, but make sure you have Python 3 and Numpy installed.

### Example command

```bash
python phasepaint.py -i input.tsv.gz -o ouptput.phasepaint_I20.tsv.gz --threads 20 --max_iterations 20 --ignore_first_n_columns 2 --n_best_scores 3
```

### Input and output formats

Both the input and output files are tab-separated text files (that can be gzipped), with the structure:

```
CHROM    POS    Ind1A   Ind1B   Ind2A   Ind2B    Ind3A   Ind3B
chr01    1001     2       1       1       2        1       3
chr01    1055     2       1       1       2        1       1
chr01    1190     2       1       1       2        1       1
chr01    1243     2       1       1       2        1       3
chr01    1407     1       2       1       2        1       3
chr01    1415     1       2       1       2        1       3
chr01    1580     1       2       1       2        1       3
```

Note:
* The first row is ignorred.
* The first two columns are ignored by default. Modify this with option `--ignore_first_n_columns`.
* the script assumes all tracts are from the same chromosome, so you should separate chromosomes into different input files.

### How it works

Looking at the example above, we can see a probable phase swith error in Individual 1 (see 4th and 5th rows down). `phasepaint` aims to find a state that maximises similarity among the haplotypes observed. In this case, we can see that this would be achieved by switching the phase in Individual 1 from position 1407 onwards. Looking next at Individual 3, we see another case where ancestry changes. However, in this case, there is no way to change the phase such that it would increase similarity among haplotypes, so Individual 3 should be left unchanged.

The above logic is applied by `phasepaint` using a heuristic algorithm. For each diploid individual, we iterate over each row in the input, and switch the phase if this would cause an *increase* in similarity of the two resulting haplotypes to one or more other haplotypes in the full dataset (considering all ancestry blocks up to and including the focal one). By default, we consider the average similarity of the **top three** most similar haplotypes in the rest of the dataset (excluding the focal one). This setting can be adjusted with option `--n_best_scores`. The process is performed left-to-right and then right-to-left. Once this phase correction procedure has been performed on all individuals, it is repeated again starting from the first individual, up to a total of `--max_iterations` iterations. `phasepaint` will only be reliable when you have a large sample size (ideally 100+). In principle, this logic can be applied for phasing genotypes directly (not only ancestry-painted SNPs or tracts). However, phasepaint has not been tested for that purpose.
