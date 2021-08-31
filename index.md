# GFLassoInference <img src="man/figures/PGInference_sticker.png" align="right" width="150px"/>

### What is GFLassoInference?

`GFLassoInference` is an `R` package for testing for a difference in means between a pair of connected components resulting from the graph fused lasso.

### How do I install the package?

To download the GFLassoInference package, use the code below.
```r
require("devtools")
devtools::install_github("yiqunchen/GFLassoInference")
library(GFLassoInference)
```

### Why do we need GFLassoInference?
Testing for equality in means between two groups is one of the most fundamental task in statistics with a numerous applications: comparing vaccine efficacy between treatment and control groups, detecting the shift in copy numbers in chromosomes, to name a few. If the groups are defined a priori, i.e., without making use of the observed data, then classical hypothesis tests will control the (selective) Type I error; that is, the probability of a false rejection, *given that* the null hypothesis is tested.

However, in practice, data analysts often find themselves testing the null hypothesis of equality in means for two groups that are *a function of the same data* used for testing. Consider the graph fused lasso, a widely popular optimization problem used to reconstruct underlying signals that are *piecewise constant* on a graph. The resulting solution can be segmented into groups known as the *connected components* &mdash; that is, the set of elements of the solution that share a common value, and are connected in the original graph. Suppose we want to test whether the means are the same across two connected components, then we need to carefully craft new procedures, because now the groups are *defined through data*, and a naive procedure like the z-test will lead to inflated selective Type I error.

To tackle this problem, Hyun et al. (2018) proposed a valid test by computing the $p$-value (henceforth referred to as $p_{\text{Hyun}}$), *conditional* on all outputs of the algorithm used to obtain the graph fused lasso solution. While their proposal controls the selective Type I error, $p_{\text{Hyun}}$ conditions on too much information than what is necessary to come up with the null hypothesis under consideration, therefore leading to extremely low power in practice. 

In this work, we propose an alternative $p$-value $p_{C_1,C_2}$ which conditions *only on the pair of connected components being tested*. The test based on the resulting $p$-value will have higher power than $p_{\text{Hyun}}$, while still guaranteeing selective Type I error control. 

As an example, consider the graph fused lasso on a grid graph, constructed by connecting each node to its four closest neighbors (up, down, left, right). This leads to the two-dimensional fused lasso problem, also known as the total-variation denoising when applied to an image (Rudin et al. 1992, Tibshirani and Taylor 2011). In the leftmost panel, we display the piecewise mean structure of the signal; in the middle panel, we see that both tests based $p_{\text{Hyun}}$ or $p_{C_1,C_2}$ (i.e., rejecting $H_0$ when the $p$-value is less than $\alpha$) control the selective Type I error, but the test based z-test $p_{\text{Naive}} = \mathbb{P}(|\nu^\top Y|\geq |\nu^\top y|)$ leads to inflated selective Type I error. In the rightmost panel, we see that the test based on $p_{C_1,C_2}$ higher power than that based on $p_{\text{Hyun}}$. More details for the simulation results can be found in Section 5.2 of our manuscript. 

![](../man/figures/combined_two_d.png)
<!-- [Figure 1: (a): The piecewise mean structure of $\beta$ according to a two-dimensional grid graph. (b): Under the null hypothesis, both $p_{\text{Hyun}}$ and $p_{C_1,C_2}$ control the selective Type I error, but the z-test $p_{\text{Naive}} = \mathbb{P}(|\nu^\top Y|\geq |\nu^\top y|)$ leads to inflated selective Type I error. (c): For a given value of the effect size ($|\nu^\top\beta|/\sigma$), $p_{C_1,C_2}$ has higher power than $p_{\text{Hyun}}$. Power for both increases as a function of the effect size.] -->

### Link to additional resources
* You can learn more about the technical details in our manuscript and in the [technical details section](https://yiqunchen.github.io/GFLassoInference/articles/technical_details.html).
* You can learn more about how to use our software in the  [tutorials section](https://yiqunchen.github.io/GFLassoInference/articles/Tutorials.html).
* Finally, code and steps to reproduce the figures in our manuscript can be found in the GitHub repo [https://github.com/yiqunchen/GFLassoInference-experiments](https://github.com/yiqunchen/GFLassoInference-experiments).

### Citation

If you use `GFLassoInference` for your analysis, please cite our manuscript:

Chen YT, Jewell SW, Witten DM. (2021+) More powerful selective inference for the graph fused lasso.

### Bug Reports / Change Requests

If you encounter a bug or would like to make a change request, please file it as an issue [here](https://github.com/yiqunchen/GFLassoInference/issues).

### References

Chen YT, Jewell SW, Witten DM. (2021+) More powerful selective inference for the graph fused lasso

Fithian W, Sun D, Taylor J. (2014) Optimal Inference After Model Selection. arXiv:1410.2597 [mathST]. 

Hyun S, G’Sell M, Tibshirani RJ. (2018) Exact post-selection inference for the generalized lasso path. Electron J Stat.

Lee J, Sun D, Sun Y, Taylor J. Exact post-selection inference, with application to the lasso. Ann Stat. 2016;44(3):907-927. doi:10.1214/15-AOS1371

Rudin L, Osher S, Fatemi E. Nonlinear total variation based noise removal algorithms. Physica D. 1992;60(1):259-268. doi:10.1016/0167-2789(92)90242-F

Tibshirani RJ, Taylor J. The solution path of the generalized lasso. Ann Stat. 2011;39(3):1335-1371. doi:10.1214/11-AOS878


