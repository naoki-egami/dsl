## dsl: Design-based Supervised Learning

------------------------------------------------------------------------

### Overview

R package `dsl` implements design-based supervised learning (DSL)
proposed in [Egami, Hinck, Stewart, and Wei
(2023)](https://naokiegami.com/paper/dsl.pdf).

DSL is a general estimation framework for using predicted variables in
statistical analyses. The package is especially useful for researchers
trying to use large language models (LLMs) to annotate a large number of
documents they analyze subsequently. DSL allows users to obtain
statistically valid estimates and standard errors, even when LLM
annotations contain arbitrary non-random prediction errors and biases.

Note: While [Egami, Hinck, Stewart, and Wei
(2023)](https://naokiegami.com/paper/dsl.pdf) focused on settings where
the outcome variable is text-based and requires prediction, this package
allows users to incorporate more general cases where any subset of the
outcome and independent variables are text-based, relying on new results
derived in our forthcoming paper.

To learn how to use the package, please start with **[Get Started
Page](http://naokiegami.com/dsl/articles/intro.html).**

### Installation Instructions

You can install the most recent development version using the `devtools`
package. First you have to install `devtools` using the following code.
Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install
`dsl`:

``` r
library(devtools)
install_github("naoki-egami/dsl", dependencies = TRUE)
```

### Information

**Authors:**

-   [Naoki Egami](https://naokiegami.com) (Maintainer)

-   [Musashi Hinck](https://muhark.github.io/about)

-   [Brandon M. Stewart](https://bstewart.scholar.princeton.edu/)

-   [Hanying Wei](https://polisci.columbia.edu/content/hanying-wei)

**Reference:**

-   [Egami, Hinck, Stewart, and Wei.
    (2023)](https://naokiegami.com/paper/dsl.pdf). Using Imperfect
    Surrogates for Downstream Inference: Design-based Supervised
    Learning for Social Science Applications of Large Language Models,
    Advances in Neural Information Processing Systems (NeurIPS).
