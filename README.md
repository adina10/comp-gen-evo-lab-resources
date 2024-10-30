# Examples of simulation code in C and figure generating in R

This repository contains:
- Simulation code in C used in my senior thesis titled "Epistasis modulates balanced polymorphism in changing environments"
- Example R code for one of my figures.
- Always feel free to ask any questions on this code or your own code.

## Descriptions

### 1. `HCS2.c`
- **Description**: This file implements Heterogeneous Cyclic Selection in a subdivided population of haploid individuals with 2 loci with epistatic interactions. It takes in a "pars1" text file for parameters and a "trajectory" file of points on the sin curve that is used to model selective pressures.

### 2. `HCS2r.c`
- **Description**: This file is an extension of HCS2.c (above) with the same input files, but implements recurrent mutation instead of single mutant introduction.

### 3. `tra.c`
- **Description**: This file generates the "trajectory" input file for HCS2.c and HCS2r.c. It takes in a text file "parst" for parameters of the sin curve (period, total cycles, smax).

### 4. `generate_figure.R`
- **Description**: This file is what I used in RStudio to generate one of my figures. It takes in a CSV file "input" with columns A) group number from 1-9 B) number of migrants in simulation C) type of epistasis D) heterozygosity levels. I use bar graphs but there are many other ways to represent your data depending on output of your code.

### 5. `figure.png`
- **Description**: This image is an example of output from generate_figure.R.
