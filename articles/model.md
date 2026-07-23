# Models Available in DSL

**Table of Contents**

- [**`lm`: Linear Regression**](#lm)  

- [**`logit`: Logistic Regression**](#logit)  

- [**`felm`: Fixed Effects Linear Regression**](#felm)  

Note: We currently cover three models (`lm`, `logit`, `felm`). In the
next update, we are planning to add more models, including
multinomial-logit regression, poisson regression, and two-stage least
squares for the instrumental variable method.

  

### **`lm`: Linear Regression**

Researchers can implement DSL linear regression by specifying
`model = "lm"`.

``` r

library(dsl)
data("data_lm") # example data 
head(data_lm)
```

    ##          Y   pred_Y         X1         X2         X3        X4         X5
    ## 1 3.337877 3.288469 0.84678708 0.10927346 0.15178520 0.4336862 0.07082542
    ## 2 1.283851 1.353864 0.59445495 0.19281988 0.77950672 0.4498675 0.30573402
    ## 3 3.618762 3.603551 0.91183977 0.97745425 0.19939829 0.3630903 0.97453388
    ## 4       NA 1.402226 0.24259237 0.03956956 0.02819825 0.8497914 0.95381503
    ## 5       NA 3.842970 0.51535594 0.78112420 0.31758914 0.5531493 0.46124773
    ## 6       NA 0.601493 0.09942167 0.28733200 0.44238525 0.5460361 0.83561538

In this example, variable `Y` requires text annotation (they take `NA`
if a given observation is not labeled by experts). Variable `pred_Y`
represents predictions for `Y`.

``` r

out_lm <- dsl(model = "lm", 
              formula = Y ~ X1 + X2 + X3 + X4 + X5,
              predicted_var = "Y",
              prediction = "pred_Y",
              data = data_lm)
```

    ## Cross-Fitting: 1/10..2/10..3/10..4/10..5/10..6/10..7/10..8/10..9/10..10/10..

| Argument | Description |
|:---|:---|
| `model` | A regression model. `dsl` currently supports `lm` (linear regression), `logit` (logistic regression), and `felm` (fixed-effects regression). |
| `formula` | A formula used in the specified regression model. |
| `predicted_var` | A vector of column names in the data that correspond to variables that need to be predicted. |
| `prediction` | A vector of column names in the data that correspond to predictions of `predicted_var`. |
| `data` | A data frame. The class should be `data.frame`. |

``` r

summary(out_lm)
```

    ## ==================
    ## DSL Specification:
    ## ==================
    ## Model:  lm
    ## Call:  Y ~ X1 + X2 + X3 + X4 + X5
    ## 
    ## Predicted Variables:  Y
    ## Prediction:  pred_Y
    ## 
    ## Number of Labeled Observations:  245
    ## Random Sampling for Labeling with Equal Probability: Yes
    ## 
    ## =============
    ## Coefficients:
    ## =============
    ##             Estimate Std. Error CI Lower CI Upper p value    
    ## (Intercept)   0.6450     0.0714   0.5050   0.7850  0.0000 ***
    ## X1            2.1293     0.0685   1.9951   2.2635  0.0000 ***
    ## X2            2.1033     0.0645   1.9769   2.2296  0.0000 ***
    ## X3            0.0448     0.0687  -0.0898   0.1794  0.5144    
    ## X4            0.0022     0.0635  -0.1223   0.1267  0.9722    
    ## X5            0.0257     0.0593  -0.0906   0.1419  0.6652    
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 95% confidence intervals (CI) are reported.

  

### **`logit`: Logistic Regression**

Researchers can implement DSL logistic regression by specifying
`model = "logit"`. It requires the same arguments as `model = "lm"`.

``` r

data("data_logit") # example data 
head(data_logit)
```

    ##    Y         X1 pred_Y    pred_X1         X2          X3          X4
    ## 1  0 -0.1443742      0  0.9486148 -0.2079232 -0.06932595  1.65748976
    ## 2  1 -1.2864079      1 -1.8766233 -1.6961906 -0.27274691  1.54678252
    ## 3  1 -1.4665386      1 -2.4866911 -0.6212224 -0.70787173 -0.27913330
    ## 4  0 -1.3554115      1  0.5138155  1.2299990 -0.86519213  0.18177847
    ## 5 NA         NA      1 -2.1116582  0.6624512  1.19504665 -0.02724476
    ## 6 NA         NA      1 -0.6391120 -0.8971608  1.20072203 -0.51241357
    ##            X5          X6         X7        X8           X9        X10
    ## 1  0.96825927  0.05389556 -2.0762690 0.2653279  1.458743948  0.4127003
    ## 2 -0.46626141  0.73182950 -1.7677551 1.0790758 -0.003172833 -0.5313550
    ## 3 -0.26973904 -1.75925076 -0.1585703 0.6823146  0.050500221 -1.2679467
    ## 4  0.41947288 -0.32375937  0.3226157 0.5131867 -2.065367939  1.6111660
    ## 5 -0.94927647  1.01808653  0.7608077 1.6830265 -1.907959205 -0.6979879
    ## 6  0.09162689  1.15086391 -0.7917356 0.9765960 -0.144924483 -1.2848539

In this example, variables `Y` and `X1` require text annotation (they
take `NA` if a given observation is not labeled by experts). Variables
`pred_Y` and `pred_X1` represent predictions for `Y` and `X1`. As in
this example, `dsl` can handle cases when multiple variables require
text annotations.

``` r

out_logit <- dsl(model = "logit", 
                 formula = Y ~ X1 + X2 + X4,
                 predicted_var = c("Y", "X1"),
                 prediction = c("pred_Y", "pred_X1"),
                 data = data_logit)
```

    ## Cross-Fitting: 1/10..2/10..3/10..4/10..5/10..6/10..7/10..8/10..9/10..10/10..

``` r

summary(out_logit)
```

    ## ==================
    ## DSL Specification:
    ## ==================
    ## Model:  logit
    ## Call:  Y ~ X1 + X2 + X4
    ## 
    ## Predicted Variables:  Y, X1
    ## Prediction:  pred_Y, pred_X1
    ## 
    ## Number of Labeled Observations:  487
    ## Random Sampling for Labeling with Equal Probability: Yes
    ## 
    ## =============
    ## Coefficients:
    ## =============
    ##             Estimate Std. Error CI Lower CI Upper p value    
    ## (Intercept)  -0.7386     0.1000  -0.9346  -0.5426  0.0000 ***
    ## X1            0.2208     0.1044   0.0162   0.4255  0.0344   *
    ## X2            0.0242     0.0974  -0.1667   0.2151  0.8038    
    ## X4            0.2742     0.1022   0.0739   0.4745  0.0073  **
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 95% confidence intervals (CI) are reported.

  

### **`felm`: Fixed Effects Linear Regression**

Researchers can implement DSL fixed effects linear regression by
specifying `model = "felm"`.

``` r

data("data_lm")  # example data 
head(data_felm)
```

    ##     state year  log_gsp pred_log_gsp log_pcap   log_pc unemp
    ## 1 ALABAMA 1970       NA     9.738913 9.617981 10.48553   4.7
    ## 2 ALABAMA 1971 10.28790    10.287899 9.648720 10.52675   5.2
    ## 3 ALABAMA 1972 10.35147    10.351469 9.678618 10.56283   4.7
    ## 4 ALABAMA 1973 10.41721    10.417209 9.705418 10.59873   3.9
    ## 5 ALABAMA 1974 10.42671    10.426706 9.726910 10.64679   5.5
    ## 6 ALABAMA 1975 10.42240    10.422400 9.759401 10.69130   7.7

In this example, variable `log_gsp` requires text annotation (it takes
`NA` if a given observation is not labeled by experts). Variable
`pred_log_gsp` represents predictions for `log_gsp`.

#### One-way Fixed Effects

To implement DSL one-way fixed effects regression, users can set
`fixed_effect = "oneway"`. Use `index` to denote which column defines
fixed effects. Users can also cluster standard errors by specifying a
variable name in argument `cluster` (this `cluster` argument is also
avaiable for other models).

``` r

out_felm_one <- dsl(model = "felm", 
                    formula = log_pcap ~ log_gsp + log_pc + unemp,
                    predicted_var =  "log_gsp",
                    prediction = "pred_log_gsp",
                    fixed_effect = "oneway", 
                    index = c("state"), 
                    cluster = "state",
                    data = data_felm)
```

    ## Cross-Fitting: 1/10..2/10..3/10..4/10..5/10..6/10..7/10..8/10..9/10..10/10..

| Argument | Description |
|:---|:---|
| `fixed_effect` | A type of fixed effects regression you run. `oneway` (one-way fixed effects) or `twoways` (two-way fixed effects). |
| `index` | A vector of column names specifying fixed effects. When `fixed_effect = oneway`, it has one element. When `fixed_effect = twoways`, it has two elements, e.g., `index = c("state", "year")`. |
| `cluster` | A column name in the data that indicates the level at which cluster standard errors are calculated. Default is `NULL`. |

``` r

summary(out_felm_one)
```

    ## ==================
    ## DSL Specification:
    ## ==================
    ## Model:  felm (oneway)
    ## Call:  log_pcap ~ log_gsp + log_pc + unemp
    ## Fixed Effects:  state
    ## 
    ## Predicted Variables:  log_gsp
    ## Prediction:  pred_log_gsp
    ## 
    ## Number of Labeled Observations:  334
    ## Random Sampling for Labeling with Equal Probability: Yes
    ## 
    ## =============
    ## Coefficients:
    ## =============
    ##         Estimate Std. Error CI Lower CI Upper p value    
    ## log_gsp   0.0031     0.0020  -0.0007   0.0069  0.1117    
    ## log_pc    0.5357     0.0419   0.4537   0.6178  0.0000 ***
    ## unemp     0.0067     0.0021   0.0025   0.0108  0.0016  **
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 95% confidence intervals (CI) are reported.
    ## Standard errors are clustered by state.

#### Two-way Fixed Effects

To implement DSL two-way fixed effects regression, users can set
`fixed_effect = "twoways"`. Use `index` (a vector of length 2) to denote
which columns define fixed effects.

``` r

out_felm_two <- dsl(model = "felm", 
                    formula = log_pcap ~ log_gsp + log_pc + unemp,
                    predicted_var =  "log_gsp",
                    prediction = "pred_log_gsp",
                    fixed_effect = "twoways", 
                    index = c("state", "year"), 
                    cluster = "state",
                    data = data_felm)
```

    ## Cross-Fitting: 1/10..2/10..3/10..4/10..5/10..6/10..7/10..8/10..9/10..10/10..

``` r

summary(out_felm_two)
```

    ## ==================
    ## DSL Specification:
    ## ==================
    ## Model:  felm (twoways)
    ## Call:  log_pcap ~ log_gsp + log_pc + unemp
    ## Fixed Effects:  state and year
    ## 
    ## Predicted Variables:  log_gsp
    ## Prediction:  pred_log_gsp
    ## 
    ## Number of Labeled Observations:  334
    ## Random Sampling for Labeling with Equal Probability: Yes
    ## 
    ## =============
    ## Coefficients:
    ## =============
    ##         Estimate Std. Error CI Lower CI Upper p value    
    ## log_gsp   0.0025     0.0017  -0.0009   0.0060  0.1445    
    ## log_pc    0.4275     0.0748   0.2809   0.5741  0.0000 ***
    ## unemp     0.0093     0.0034   0.0027   0.0159  0.0055  **
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 95% confidence intervals (CI) are reported.
    ## Standard errors are clustered by state.
