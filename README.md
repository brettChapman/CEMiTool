# CEMiTool [![Build Status](https://travis-ci.org/csbl-usp/CEMiTool.svg?branch=master)](https://travis-ci.org/csbl-usp/CEMiTool)[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

A modified version of the Co-Expression Module Identification Tool (CEMiTool) (modified by Brett Chapman)

### Docker Image
For more information about using the official CEMiTool Docker image, see [here](docker/example.md)
The Dockerfile has been modified by Brett Chapman to utilise this GitHub repo

### Updates to this version of CEMiTool

New parameters have been added to the CEMiTool.R executable in exec/ and the script has been significantly updated with added functionality.

In addition the Beta value can now be set manually or forced if R squared is below 0.8, as explained in point 6 here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html.

Usage:
```--top_hubs N (default: 10)``` to output a list of the top N hub genes per module, and a filtered expression matrix is also output for use in other tools if desired.
```--top_hubs_interact N (default: 10)``` to use the top N hub genes as the basis of the interaction network.
```--set-beta N (default: 'None')``` set the Beta value manually.
```--force-beta N (default: 'None')```  force the Beta value if R squared is found to be below 0.8.

Another parameter included which was left out from the original executable is ```--cor-function=<corfunc>```

A python script has also been included which converts GFF to a GMT file if wanting to generate a gene set for ORA analysis. Modification of the python script may be nessary depending on the formatting of your GFF file. A test.gff file has been included in tests/.

Please use this version of CEMiTool if you would like to run CEMiTool on a cluster with the executable through Docker or Singularity instead of RStudio, and if you would like to see a list of the top hub genes, use these genes for the interaction network (if you don't have a interactome of your own), and would like to slice the original expression matrix by hub genes and transformed using VST.
