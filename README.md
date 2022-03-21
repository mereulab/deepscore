![alt text](https://github.com/pabloswfly/DeepSCore/blob/master/logo_thin.png?raw=true)

# DeepSCore
### A multi-language multi-omics deep learning model for automatic single-cell label transfer.

- **/tutorials** directory includes jupyer notebook showing how to prepare and run DeepScore in both R and Python.
- **/utils** directory includes convenient extra functions that facilitates working with different fortmats of single-cell data.
- **/R** directory includes the source code to run DeepScore in R.
- **/python** directory includes the source code to run DeepScore in Python.

Even though the package is fully functional, it is still in development and it has not been packaged and released in a public repository. 
To implement deepScore in your analyses, run:

## R

```R
setwd("<dirpath>/deepscore/R")
source("deepscore.R")
source("marker_analysis.R")
```

In R, the following dependencies are needed:
- seurat / signac
- reshape2
- dplyr
- ggplot2
- tensorflow / keras

## Python

```Python
import sys
sys.path.append("<dirpath>/deepscore/python")
from deepscore import DeepScore
from marker_analysis import *
```

In R, the following dependencies are needed:
- scanpy / episcanpy
- anndata
- seaborn
- matplotlib
- pandas
- numpy
- tensorflow / keras
