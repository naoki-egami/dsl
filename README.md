# dsl: Design-based Supervised Learning

We explain the details of the package on the package website
(<http://naokiegami.com/dsl>).

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
