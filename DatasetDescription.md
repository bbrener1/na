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

## Heatmaps for fun and profit

Ok, so we've looked at 20 histograms of gene expression. Perhaps we were unlucky and got a weird selection? What are some other meaningful ways of looking at boatloads of organized data to try to find if there are similar ways in which it can behave? This sounds like a job for heatmaps and clustering!

If we simply wanted to look at the totality of the matrix without having to squint really hard at tiny numbers, a simple way to do so would be to color-code each value in a logical way, shrink it all down and look at the overall color pattern. Even by itself this can be informative to some degree. Below is a heat map, each row of pixels is a single cell, and each column is a single gene, so pixel[x,y] is the color of cell x expression of gene y (in whatever order they appear in the header of our file).

![](https://github.com/bbrener1/na/blob/master/figures/raw_expression_heatmap.png "An unorganzied heatmap of gene expression")

Ok, so far we sort of see some of the things we already expected recaptured, some genes have highish overall expression across all cells, some genes have barely any expression anywhere. This is the vertical lines. No obvious horizontal lines yet, but variability across cells isn't as great as between genes. At the moment though, the ordering along both axes is pretty much random.

Can we put these axies in a useful or interesting order? If we are assuming that there should be patterns to the gene expression that recurr among different cells, a typical approach to finding that pattern would be something like heirarchal clustering, which fortunately SciPy can do for us. 

![](https://github.com/bbrener1/na/blob/master/figures/doubly_clustered_raw_genes.png "Heirarchically clustering genes and cells")

Now let's take a look at what we got from heirarchal clustering. First off we should note that the dendrograms look a little weird. Heirarchal clustering is heavily dependent on the distance metrics and centering procedures for its final output. Now unfortunately, the most conventional distance metric, namely eucledian distance, was giving me trouble because it was reaching into too much recursive depth and breaking Scipy. But that would imply to us that it was doing something weird and dumb anyway, so for the moment let's leave that aside. Empirically, the most interesting results come around when use we  the "average" cluster center method, cosine distances for cells, and correlation distances for genes. That's the heatmap you see. 

These types of heatmaps are also fairly conventional. The squares you see are blocks of genes that are up or down regulated in specific cell types, ie a horizontal band is a cell type, and vertical bands are the specific patterns of expression of gene blocks in those cell types. Now notice something more interesting. If you remember our earlier histograms, we saw that gene expression across the whole cell population was a bimodal distribution, but it looked normalish at those two peaks. But note that in this diagram, we don't necessarily see that much "noise". There is a little speckling for highly expressed genes, but the TV-static pattern of the unordered heatmap is gone. 

On the one hand, we shouldn't be surprised that clustering brings together similar cells, we're still doing well and our data makes sense.

On the other hand, it suggests to us a slightly less obvious point. When we see a random-looking distribution in gene expression values relative to some factor (such as gene vs gene scatter plot), we might be tempted to think of the fact that distribution looks like a bell curve or something similar as being the result of technical noise or "just biology being biology". On the other hand, looking at how well ordered this heatmap looks, perhaps the degree of randomness in the average genetic expression bell curve is smaller than it appears?

## Cell Types?

Now, it might be useful for us to save this sorted version of the plot as well as the re-ordered cell and gene labels. 

## Deviation matrices

Consider local pseudo-cells

## Gene-gene relationships/Prediction

I'm going to skip ahead a little bit to investigate a question that I am currently working on, which is the ability of genes to predict each other's behavior. 

We have many gene-gene comparisons to make when we are working with 4773 genes (specifically 22 million comparisons) to determine whether the expression of one gene influences the expression of the other. Determining whether this is the case is a non-trivial thing, and one of the first steps we should take is, again, just looking at some data. 

At a most basic level we are sort of hoping that a linear, exponential, or polynomial relationship exists between the expression of one gene and the expression of another gene. We can also consider relationships that have a more complicated nature, but generally, the more complex a relationship between two genes, the worse our odds are of inferring it accurately. So let's ask some basic questions first. 

First off, when we were looking at the potential expression values that genes exhibited, we looked at a set of 20 histograms to see what kinds of gene expression values genes could take on. Let's do something similar now, but generate many random plots of gene-gene scatterplots. Considering that gene expression is bimodal, but the average genes should bear relatively little influence on each other, generally we should see four blobs on our scatter plots, one in each corner of the plot, of relatively even intensities. In genes that are able to predict each others expression we should see uneven sizes between the blobs, but it's unclear how often we should encounter such relationships. (Actually we'll answer that question a little later, but for now let's just do something quick and dirty)

![](https://github.com/bbrener1/na/blob/master/figures/gene_scatter_gigaplex.png "Various randomly chosen genes scattered against each other")

Ok, so far not so interesting, most genes don't predict each other's behavior very well. When we see non-diagonal patterns like the ones that are common on this plot, mainly it just means one of the genes is not very bimodal. (Interesting to see that Ttc30b is tri-modal, this is kind of unexpected. If this figure has been changed since then check it out sometime, it's pretty weird.) 

But there are obviously genes that correlate to each other very well, as we saw from our clustered heatmap. How do we get to see what their behavior looks like? Well the clustering parameter for genes on the heatmap was correlation, so first we'd like to know which genes correlate to which, and what kinds of correlations we see in general. 

Numpy provides us with the "corrcoef" function which does this quickly and helpfully by computing the Pearson correlation of each row in a matrix to each other row in the matrix. 

At face value, we see that low or negligible correlations are present between most genes (note that the plot is logarithmic), but some genes do have substantial correlation. 

![](https://github.com/bbrener1/na/blob/master/figures/correlation_histogram.png "Histogram of all possible correlation pairings")

Let's see what a random picking of correlated genes looks like when scattered against each other first, and then we can look at the correlations between a slice of the clustered heatmap:

![](https://github.com/bbrener1/na/blob/master/figures/correlated_scatter_gigaplex.png "Scatterings of correlated genes")

Ok, now we are seeing some more interesting plots. There are about 34000 gene-gene pairings that correlate to each other at a rate higher than .5, out of 22 million total possible correlations. 

On the other hand, if we consider correlations between .1 and .5 (still probably pretty significant, but much more common), we find that there are 3.2 million such correlations (nearly 10% of all possible gene-gene pairings), and they look something like this:

![](https://github.com/bbrener1/na/blob/master/figures/correlated_scatter_gigaplex(intermediate).png "Scatterings of (weakly) correlated genes")

Ok, so what do we see? We are observing a still bimodal sort of interaction, but essentially you see certain genes in concurrent states more often than not. Finding really unambiguous relationships though, isn't very common.

This serves as encouraging news, parsing through a network of 3.2 million gene-gene interactions would be difficult, and 34k gene-gene interactions sounds like a plausible amount if we consider scale-free network topologies as being a good representation of human genetic regulation. Unfortunately if we actually look at the correlations observed, it becomes clear that most of them are probably not direct genetic regulatory relationships but merely what they sound like, correlations. 

Let's examine this in a bit more detail by looking at some of the top of the gene-gene pairings that have the top correlations: 

> Rhd	Slc38a5	0.877749294748  
> Atp1b2	Tspo2	0.879771717678  
> Tspo2	Atp1b2	0.879771717678  
> Birc5	Ccna2	0.881493438329  
> Ccna2	Birc5	0.881493438329  
> Lars2	Gm15564	0.881922602514  
> Gm15564	Lars2	0.881922602514  
> Hba-a1	Hba-a2	0.883262737553  
> Hba-a2	Hba-a1	0.883262737553  
> Igkv14-126	Ighv11-2	0.884403538866  
> Ighv11-2	Igkv14-126	0.884403538866  
> Ctsg	Mpo	0.885019347892  
> Mpo	Ctsg	0.885019347892  
> Rhd	Sphk1	0.885555721254  
> Sphk1	Rhd	0.885555721254  
> Rhag	Rhd	0.886360349907  
> Rhd	Rhag	0.886360349907  
> Ltf	Retnlg	0.888045084374  
> Retnlg	Ltf	0.888045084374  
> Car1	Aqp1	0.888748995754  
> Aqp1	Car1	0.888748995754  
> Car1	Tspo2	0.889320897323  
> Tspo2	Car1	0.889320897323  
> Cldn13	Rhd	0.890933003973  
> Rhd	Cldn13	0.890933003973  
> Atp1b2	Car1	0.891350787339  
> Car1	Atp1b2	0.891350787339  
> Klf1	Aqp1	0.891865584038  
> Aqp1	Klf1	0.891865584038  
> Car1	Klf1	0.896917892665  
> Klf1	Car1	0.896917892665  
> Jchain	Igkc	0.909690216796  
> Igkc	Jchain	0.909690216796  
> Atp5g1	Gm10039	0.943769882157  
> Gm10039	Atp5g1	0.943769882157  
> Gm13461	Ran	0.947917392347  
> Ran	Gm13461	0.947917392347  
> Hsp90aa1	Gm12346	0.950589822034  
> Gm12346	Hsp90aa1	0.950589822034  
> Cd63	Cd63-ps	0.969345523311  
> Cd63-ps	Cd63	0.969345523311  
> Igkv4-50	Ighv9-1	0.990570856454  
> Ighv9-1	Igkv4-50	0.990570856454  

Ok, so top correlations are two immunoglobulin subunits, which is great news because these should rarely have non-1-to-1 stoichiometry. Cd63 is correlated to its pseudogene. At position 3 we are already getting into something interesting. GM12346 is a predicted gene of unknown function that is apparently tightly regulated by Hsp90. Position 4 is again a gene of unknown function, though this time annotated as a probable reverse-transcribed pseudogene, tightly correlated with Ran. Interesting. Brief investigation implies that this pseudogene is located in the middle of LRP1B, a low density lipoprotein that serves as a tumor supressor (???), and is a common target for papilomaviridae, which frequently induces cancer. This mouse lucked out with the excellent healthcare it got. Ran interaction is intersting though. Ran is a decent marker for cell cycle, so perhaps this pseudogene is cell-cycle correlated? Spot 5, more pseudogenes. This presudogene appears to be a reverse-transcribed mRNA of an ATP synthase component, and it correlates with the expression of the actual ATP synthase componenet. This may be an alignment artifact, actually. Spot 6 is more immunoglobulin componenets correlating to each other.

Now, however, we are getting into the good stuff. Here we see that KLF1 corresponds well to a couple of things that RBCs really need. Both Aquaporin and Car1 are important componenets of RBCs and this is probably the first genuine regulatory relationship we have seen so far. This is the type of result that we are really hoping to see. But, out of curiosity, let's try to take a look at the correlation between Car1 and Aqp1. Oh... It's like 10 lines up the list, and it's .88. This is the type of problem that we're going to run into a lot. If a really meaty regulatory factor like Klf1 has a tight regulatory relationship with two things, those two things will look like they're acting on each other as well. 

So what's a girl to do? This is one of the central questions of gene regulation studies. Many people throw up their hands, call these types of triplets (or more than triplets) "gene blocks", throw them up in databases, and let biologists sort it all out. To be fair, manually looking at some of these blocks, it can sometimes be obvious which transcription factor is really controlling a group of proteins. Other times, these relationships get less obvious and more muddled. 

In general, you will observe something like this: 
![](https://github.com/bbrener1/na/blob/master/figures/gene_covariance_simple.png "Heatmap of gene-gene correlations")


We can examine some of these relationships by looking at the checker patterns that jump out at us from the clustering procedure. If we look more closely at the specific checks, namely what genes and cells go into them, we can examine some of the obvious correlations between these genes and the cells they are active in, eg:

Another thing that we can do is to attempt to partition the cells in the same way. It stands to reason that if blocks of genes that co-express represent specific cellular states, then cells that occupy those states should have some degree of internal organization also. This is a cell covariance matrix:

![](https://github.com/bbrener1/na/blob/master/figures/gene_covariance_simple.png "Heatmap of cell-cell correlations")

How much information can we attribute to a simple analysis of where cells rest in this space?
A fair amount, KNN imputation of different cell-gene values generally has a correlation of ~.6, although the mean squared error is quite high at ~6, which means on average we are something like 4 orders of magnitude off about the expression value of any given gene. 

Let's take a look at how KNN and some other imputation methods perform at computing cell-gene values that we hide from them.



Other approaches try to look at things like mutual information. ARACNE is a good example, esssentially, by examining the triplets like this and computing mutual information three ways, they purport to discover only the most direct relationships, since theoretically the mutual information between two genes will be greatest when there is a direct relationship between them, and will be less if there is a mediating gene. 



I had a slightly different train of thought. In order to look at that, let's take alook at some of those cells we so helpfully grouped together earlier. 

