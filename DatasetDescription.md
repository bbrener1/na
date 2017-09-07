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

