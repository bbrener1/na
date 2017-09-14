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

Another interesting way to visualize this type of data is to try to use both of these pieces of information simultaneously.

In a sense, do we generally see higher expression values for genes in general when size of the recovered transcriptome for the cell is larger? Or no? A way to examine this question is to look at a plot of gene expression level plotted against the size of the transcriptome for the cell the gene came from. Here we get something like this:

![](https://github.com/bbrener1/na/blob/master/figures/trans_size_vs_exp.png "A Weird Graph")

I'm not entirely sure what to make of this data yet, but it's not obviously a diagonal so that's mixed news. Generally it looks like getting a larger transcriptome is a consequence either of only very few very high expression genes (wouldn't be apparent on this plot since I made most dots mostly transparent) or capturing more rare genes at higher levels. I'll think about this problem and come back to it later.

## Ok so how do we usefully talk about 4k genes at once? 

Well, we have a broad description of the the dataset, at least we are getting a feel for some of the numbers. 

First things first, let's try to get a feel for what the distributions of expression values for these genes could look like. There are a few ways of doing this. Step 1 might to just grab 20 histograms at random. Let's try that. 

![](https://github.com/bbrener1/na/blob/master/figures/gene_histogram_gigaplex.png "Many gene expression histograms")
(Note, this particular figure isn't deterministic to the code that generated it, but if you go back and generate histograms for the genes that do appear here, they should look identical)

Pardon the janky formatting, lower axies are gene expression values, y axies are frequencies in log space. We're not super interested in the particular details anyway, so don't look too hard. The thing to notice is that generally speaking gene expression looks bimodal. This isn't an entirely unexpected finding, there are actually many cell collections where it's the case that the average gene will be "high", "low" or "off" , but it's good to confirm and gives us a bit of a concept of what to expect.

Let's hang on to this idea and we'll examine it and its implications later.

TODO: Cluster uniform histograms also, see what happens

Ok, so we've looked at 20 histograms of gene expression. Perhaps we were unlucky and got a weird selection? What are some other meaningful ways of looking at boatloads of organized data to try to find if there are similar ways in which it can behave? This sounds like a job for heatmaps and clustering!

If we simply wanted to look at the totality of the matrix without having to squint really hard at tiny numbers, a simple way to do so would be to color-code each value in a logical way, shrink it all down and look at the overall color pattern. Even by itself this can be informative to some degree. Below is a heat map, each row of pixels is a single cell, and each column is a single gene, so pixel[x,y] is the color of cell x expression of gene y (in whatever order they appear in the header of our file).

![](https://github.com/bbrener1/na/blob/master/figures/raw_expression_heatmap.png "An unorganzied heatmap of gene expression")

Ok, so far we sort of see some of the things we already expected recaptured, some genes have highish overall expression across all cells, some genes have barely any expression anywhere. This is the vertical lines. No obvious horizontal lines yet, but variability across cells isn't as great as between genes. At the moment though, the ordering along both axes is pretty much random.

Can we put these axies in a useful or interesting order? If we are assuming that there should be patterns to the gene expression that recurr among different cells, a typical approach to finding that pattern would be something like heirarchal clustering, which fortunately SciPy can do for us. 

![](https://github.com/bbrener1/na/blob/master/figures/doubly_clustered_raw_genes.png "Heirarchically clustering genes and cells")

Now let's take a look at what we got from heirarchal clustering. 
