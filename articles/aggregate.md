# What if the Unit of Analysis is not the same as the Unit of Labeling?

  

**On this page, we explain how to setup data for `dsl` when the unit of
analysis is not the same as the unit of labeling.**

### **Overview**

In some applications, the unit of analysis is not the same as the unit
of labeling. Researchers often annotate each document, but the unit of
analysis might be some aggregates of documents. For example, users might
code whether each online post by political candidates mentions an
economic policy: the unit of labeling is at the online post level. But
researchers might be interested in how the proportion of posts
mentioning an economic policy varies between different candidates: the
unit of analysis is at the candidate level. Here, each candidate has
multiple posts, and the main text-based variable is defined as the
proportion of online posts mentioning an economic policy for each
candidate.

In such applications, how can users apply DSL? In general, DSL should be
applied at the unit of analysis (i.e., the candidate-level data). Then,
how can we prepare the candidate-level data from the post-level data
that are expert-coded? We need two simple steps: (1) aggregating
expert-coding and (2) aggregating sampling probabilities

  

#### **1. Aggregating Expert-Coding**

First, we discuss how to aggregate expert coding at the post-level data
up to the candidate-level data.

We start with the post level data.

``` r

library(dsl)
data("data_post") # example data 
head(data_post, n = 7)
```

    ##   cand_id post_id economy pred_economy post_sample_prob
    ## 1       1     1_1       0            0              0.3
    ## 2       1     1_2       0            1              0.3
    ## 6       2     2_1      NA            0              0.3
    ## 7       2     2_2      NA            1              0.3
    ## 3       3     3_1       1            1              0.3
    ## 4       3     3_2       0            0              0.3
    ## 5       3     3_3      NA            0              0.3

In this data, variable `economy` represents expert-coded labels about
whether a given post mentions an economic policy. The first candidate
has two posts that are all labeled, and for this candidate, we can
compute the proportion of online posts mentioning an economic policy to
be `0`. For the second candidates, none of the posts are labeled, so
this second candidate is not labeled. For the third candidate, two of
the posts are labeled. As long as posts are randomly selected for expert
annotations, researchers can use the mean of labeled posts to estimate
the proportion of online posts mentioning an economic policy to be `0.5`
for this particular candidate.

Variable `pred_economy` represents text labels predicted by a
user-specified automated text annotation method. Variable
`post_sample_prob` is the post-level sampling probability for labeling.
In this example, all the posts have the equal probability of being
sampled with 30%.

In general, when posts are randomly selected for expert annotations,
users can compute the expert-coded `mean_economy` variable at the
candidate-level by ignoring `NAs` and computing the proportion of online
posts mentioning an economic policy just focusing on labeled documents.
For candidates that do not have any labeled posts, those candidates will
be recorded as non-labeled in the candidate-level data. For variable
`pred_economy`, users can directly compute the mean within each
candidate.

Users can implement this step by simply using
[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) and
[`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html)
functions in `tidyverse` or using
[`tapply()`](https://rdrr.io/r/base/tapply.html).

``` r

library(tidyverse)
economy_cand <- data_post %>% 
  group_by(cand_id) %>%
  summarize(mean_economy = mean(economy, na.rm = T),
            mean_pred_economy = mean(pred_economy, na.rm = T))

# Turn NaN to NA (get ready for `dsl`)
economy_cand$mean_economy[is.nan(economy_cand$mean_economy)] <- NA 

head(economy_cand)
```

    ## # A tibble: 6 × 3
    ##   cand_id mean_economy mean_pred_economy
    ##     <dbl>        <dbl>             <dbl>
    ## 1       1          0               0.5  
    ## 2       2         NA               0.5  
    ## 3       3          0.5             0.333
    ## 4       4          0               0.5  
    ## 5       5          0               0.25 
    ## 6       6         NA               0.667

At the candidate-level, variable `mean_economy` is the expert-coded
proportion of posts mentioning an economic policy, and variable
`mean_pred_economy` is the predicted proportion of posts mentioning an
economic policy.

Here, we are using the mean to aggregate posts within each candidate,
but users can use any function of their choice to define their
text-based variable. In general, as long as posts are randomly selected
for expert annotations, users simply need to apply their definition of
text-based variables for labeled documents (ignoring non-labeled
documents) and treat candidates that have no labeled posts as
non-labeled.

  

#### **2. Aggregating Sampling Probabilities**

Users also need to translate the sampling probability from the
post-level to the candidate-level.

In particular, the sampling probability of having expert-coding for each
candidate is the probability of having expert-coding on at least one
document within each candidate. When users randomly sample documents,
they need to take into account different numbers of posts within each
candidate. This is because if a given candidate has more posts, she is
more likely to have at least one expert-coded document.

To compute this properly, users can rely on the following simple
function. Mathematically this computes the probability of having at
least one expert annotation for each candidate.

``` r

cand_sample_prob <- data_post %>% 
  group_by(cand_id) %>%
  summarize(cand_sample_prob = 1 - prod(1 - post_sample_prob))

head(cand_sample_prob)
```

    ## # A tibble: 6 × 2
    ##   cand_id cand_sample_prob
    ##     <dbl>            <dbl>
    ## 1       1            0.51 
    ## 2       2            0.51 
    ## 3       3            0.657
    ## 4       4            0.760
    ## 5       5            0.760
    ## 6       6            0.657

If users already know that their main statistical analyses are at the
candidate level, they can consider two-stage sampling, i.e., randomly
sample candidates first and then randomly sample posts within each
candidate. Using this two-stage sampling will make the calculation of
the sampling probabilities at the candidate level even simpler. They
have to just use the first stage (the probability of sampling each
candidate for expert annotations).

  

#### **3. Merging Other Candidate-level data**

To prepare the final candidate-level data, we can combine two data sets
above as well as other candidate-level data.

``` r

# First we merge (mean_economy, mean_pred_economy) and cand_sample_prob
cand_main <- merge(economy_cand, cand_sample_prob, by = "cand_id")
```

Then, users might merge some candidate-level data.

``` r

data("data_cand") # example candidate-level data 
data_cand_merged <- merge(cand_main, data_cand, by = "cand_id")
head(data_cand_merged)
```

    ##   cand_id mean_economy mean_pred_economy cand_sample_prob gender edu ideology
    ## 1       1          0.0         0.5000000           0.5100      0   2        5
    ## 2       2           NA         0.5000000           0.5100      0   4        3
    ## 3       3          0.5         0.3333333           0.6570      1   4        2
    ## 4       4          0.0         0.5000000           0.7599      1   4        3
    ## 5       5          0.0         0.2500000           0.7599      0   5        4
    ## 6       6           NA         0.6666667           0.6570      1   4        3

  

#### **4. Fitting DSL with the Candidate-level Data**

Finally, run `dsl` on the candidate-level data by specifying the
sampling probability.

``` r

out_agg <- dsl(model = "lm", 
               formula = mean_economy ~ gender + edu + ideology,
               predicted_var =  "mean_economy",
               prediction = "mean_pred_economy",
               sample_prob = "cand_sample_prob",
               data = data_cand_merged)
```

    ## Cross-Fitting: 1/10..2/10..3/10..4/10..5/10..6/10..7/10..8/10..9/10..10/10..

``` r

summary(out_agg)
```

    ## ==================
    ## DSL Specification:
    ## ==================
    ## Model:  lm
    ## Call:  mean_economy ~ gender + edu + ideology
    ## 
    ## Predicted Variables:  mean_economy
    ## Prediction:  mean_pred_economy
    ## 
    ## Number of Labeled Observations:  674
    ## Random Sampling for Labeling with Equal Probability:  No
    ## (Sampling probabilities are defined in `sample_prob`)
    ## 
    ## =============
    ## Coefficients:
    ## =============
    ##             Estimate Std. Error CI Lower CI Upper p value    
    ## (Intercept)   0.4057     0.0696   0.2693   0.5421  0.0000 ***
    ## gender       -0.0099     0.0348  -0.0780   0.0583  0.7765    
    ## edu           0.0005     0.0128  -0.0246   0.0256  0.9703    
    ## ideology     -0.0091     0.0134  -0.0353   0.0172  0.4985    
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 95% confidence intervals (CI) are reported.
