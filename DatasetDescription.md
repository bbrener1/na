# Dataset Description 

This software is intended to function with a Single-Cell RNA Sequencing dataset. The development was performed with the Vision Dataset, available here:

http://blood.stemcells.cam.ac.uk./single_cell_atlas.html

The cell/gene matrix provided by this page is broadly-speaking Log2-transformed, but uses a proprietary normalization procedure based on clustering related cells and ERCC Spikeins

The source paper for the dataset is:
https://doi.org/10.1182/blood-2016-05-716480
> Nestorowa, Sonia, Fiona K. Hamey, Blanca Pijuan Sala, Evangelia Diamanti, Mairi Shepherd, Elisa Laurenti, Nicola K. Wilson, David G. Kent, and Berthold Göttgens. "A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation." Blood 128, no. 8 (2016): e20-e31.

And the normalization procedure used to generate the normalized DE matrix provided is:
https://doi.org/10.1186/s13059-016-0947-7
> Lun, Aaron TL, Karsten Bach, and John C. Marioni. "Pooling across cells to normalize single-cell RNA sequencing data with many zero counts." Genome biology 17, no. 1 (2016): 75.	

## What's In This Dataset?

In its raw form the dataset is a 1656x4773 matrix of numeric values, each row representing cells and each column representing the (Normalized, I think?) gene counts in Log2 space. Additionally, the matrix contains fluorescence values for various FACS readouts, mostly antibodies against membrane proteins, but also including DAPI. Finally the first 3 columns are a set of 3 coordinates mapping each cell to a location in a 3D diffusion map. 

## How do the values look?

First let's take a very top-level look, how many transcripts are we seeing per cell, at least roughly?

![](https://github.com/bbrener1/na/blob/master/figures/transcript_totals_hist.png "Transcript Total Histogram")
(to obtain this figure I took the counts matrix, exponentiated each value by 2 (ie 2^counts[i,j]), and summed across each cell, then plotted resulting values in a histogram)


Ok, so the median cell yielded about 400k transcripts, and we don't really see many cells below 200k transcripts, which makes sense because most such cells were filtered through quality controls before normalization. (Presumably the ones we do see below 200k were normalized to there.

What does this look like in terms of individual gene expression?

![](https://github.com/bbrener1/na/blob/master/figures/counts_frequency.png "Gene Expression Histogram")

This is a histogram of all values present in the matrix, so you are looking at the frequency of seeing individual genes expressed at various log2 expression values. Nothing especially mind-blowing, median gene expression is approx 250 copies, presumably mostly encountered in 400k transcript cells. Perhaps about half of the values observed are below 250 copies, with many genes detected at the quite low/unreliable level of 1-2 copies. Observe the large number of raw 0 values. Of course in a log matrix a value of 0 means 1 copy, but we can probably safely round most of these values down to an actual 0. It's a common issue with single-cell RNA seq to encounter many genes that are spuriously marked as possessing 0 expression. We'll think about how to deal with this issue later though.

No pretty picture available, but roughly 40% of the cells in the matrix are zeroes. 
