# PGInference <img src="man/figures/PGInference_sticker.png" align="right" width="150px"/>

### What is PGInference?

`PGInference` is an `R` package for performing Powerful Graph fused lasso Inference.

### How do I install the package?

To download the PGInference package, use the code below.
``` r
require("devtools")
devtools::install_github("yiqunchen/PGInference")
library(PGInference)
```

### Why do we need PGInference?
Double-dipping &mdash; or more formally, generating a hypothesis based on your data, and then testing the hypothesis on that same data &mdash; renders classical hypothesis tests (e.g., z-test or Wilcoxon-rank for a difference-in-means; and in general any standard hypothesis tests) invalid, in the sense that Type I error is not controlled.

### Link to additional resources
* You can learn more about the technical details in our manuscript and in the [technical details section](https://yiqunchen.github.io/PGInference/articles/technical_details.html).
* You can learn more about how to use our software in the  [tutorials section](https://yiqunchen.github.io/SpikeInference/articles/Tutorials.html).
* Finally, code and steps to reproduce the figures in our manuscript can be found in the GitHub repo [https://github.com/yiqunchen/SpikeInference-experiments](https://github.com/yiqunchen/SpikeInference-experiments).

### Citation

If you use `PGInference` for your analysis, please cite our manuscript:

Chen YT, Jewell SW, Witten DM. (2021+) More powerful selective inference for the graph fused lasso.

### Bug Reports / Change Requests

If you encounter a bug or would like to make a change request, please file it as an issue [here](https://github.com/yiqunchen/PGInference/issues).

### References

Chen YT, Jewell SW, Witten DM. (2021+) More powerful selective inference for the graph fused lasso

Fithian W, Sun D, Taylor J. (2014) Optimal Inference After Model Selection. arXiv:1410.2597 [mathST]. 

Hyun S, G’Sell M, Tibshirani RJ. (2018) Exact post-selection inference for the generalized lasso path. Electron J Stat.

Lee JD, Sun DL, Sun Y, Taylor JE. Exact post-selection inference, with application to the lasso. Ann Stat. 2016;44(3):907-927. doi:10.1214/15-AOS1371


