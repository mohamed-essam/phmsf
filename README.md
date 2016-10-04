# phmsf
Parallel Heuristic for Minimum Spanning Forests.

##Examples

USC SIPI 2.2.13
![USC-SIPI 2.2.13 Stockton](http://i.imgur.com/eMoaxRq.png)
USC SIPI testpat.1k
![USC SIPI testpat1.k](http://i.imgur.com/Cr2jXAf.png)

##Usage

segmentImage(imageData, height, width, minimumEdgeWeight, maximumEdgeWeight, minimumRegionSize)

Parameter|Description
---------|----------
imageData|unsigned char array containing rgb values of image in a sequence {r1,g1,b1,r2,g2,b2,.....}
height, width|Height and width of image
minimumEdgeWeight|The edge weight where all edges lighter than it are automatically merged
maximumEdgeWeight|The maxmimum edge weight at which no regions are merged
minimumRegionSize|The minimum region size

