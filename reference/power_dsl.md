# Power Analysis for DSL Regression

Power Analysis for DSL Regression

## Usage

``` r
power_dsl(
  labeled_size = NULL,
  dsl_out = NULL,
  model = "lm",
  formula,
  predicted_var,
  prediction = NULL,
  data,
  cluster = NULL,
  labeled = NULL,
  sample_prob = NULL,
  index = NULL,
  fixed_effect = "oneway",
  sl_method = "grf",
  feature = NULL,
  family = "gaussian",
  cross_fit = 5,
  sample_split = 10,
  seed = 1234
)
```

## Arguments

- labeled_size:

  A vector indicating the number of labeled documents for which the
  function predicts standard errors.

- dsl_out:

  An output from function `dsl`. When this is supplied, the remaining
  arguments are overwritten by arguments specified in the output of
  `dsl`. When this is `NULL`, the function will use arguments specified
  below.

- model:

  A regression model `dsl` currently supports `lm` (linear regression),
  `logit` (logistic regression), and `felm` (fixed-effects regression).

- formula:

  A formula used in the specified regression model.

- predicted_var:

  A vector of column names in the data that correspond to variables that
  need to be predicted.

- prediction:

  A vector of column names in the data that correspond to predictions of
  `predicted_var`.

- data:

  A data frame. The class should be `data.frame`.

- cluster:

  A column name in the data that indicates the level at which cluster
  standard errors are calculated. Default is `NULL`.

- labeled:

  (Optional) A column name in the data that indicates which observation
  is labeled. It should be a vector of 1 (labeled) and 0 (non-labeled).
  When `NULL`, the function assumes that observations that have `NA` in
  `predicted_var` are non-labeled and other observations are labeled.

- sample_prob:

  (Optional) A column name in the data that correspond to the sampling
  probability for labeling a particular observation. When `NULL`, the
  function assumes random sampling with equal probabilities.

- index:

  (Used when `model = "felm"`) A vector of column names specifying fixed
  effects. When `fixed_effect = oneway`, it has one element. When
  `fixed_effect = twoways`, it has two elements, e.g.,
  `index = c("state", "year")`.

- fixed_effect:

  (Used when `model = "felm"`) A type of fixed effects regression you
  run. `oneway` (one-way fixed effects) or `twoways` (two-way fixed
  effects).

- sl_method:

  A name of a supervised machine learning model used internally to
  predict `predicted_var` by fine-tuning `prediction` or using
  predictors (specified in `feature`) when `prediction = NULL`. Users
  can run
  [`available_method()`](http://naokiegami.com/dsl/reference/available_method.md)
  to see available supervised machine learning methods. Default is `grf`
  (generalized random forest).

- feature:

  A vector of column names in the data that correspond to predictors
  used to fit a supervised machine learning (specified in `sl_method`).

- family:

  (Used when making predictions) A variable type of `predicted_var`.
  Default is `gaussian`.

- cross_fit:

  The fold of cross-fitting. Default is `5`.

- sample_split:

  The number of sampling-splitting. Default is `10`.

- seed:

  Numeric `seed` used internally. Default is `1234`.

## Value

`dsl` returns an object of `dsl` class.

- `predicted_se`: Predicted standard errors for coefficients. The first
  row shows the current standard errors for coefficients. The remaining
  rows show predicted standard errors.

- `labeled_size`: A vector indicating the number of labeled documents
  for which the function predicts standard errors.

- `dsl_out`: An output from function `dsl`.
