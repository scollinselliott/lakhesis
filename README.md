
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="right" width="150px"/> lakhesis: Consensus Seriation for Binary Data

<!-- badges: start -->
<!-- badges: end -->

The `R` package `lakhesis` provides an interactive platform and critical
measures for seriating binary data matrices through the exploration,
selection, and consensus of partially seriated sequences.

In a word, seriation (sequencing, ordination) involves putting a set of
things in an optimal order. In archaeology, seriation can be used to
establish a chronological order of contexts and find-types on the basis
of their similarity, i.e, that things come into and go out of fashion
with a peak moment of popularity (Ihm 2005). In ecology, the
distribution of a species may occur according to a preferred
environmental condition that diminishes as that environment changes (ter
Braak and Looman 1986). There are a number of R functions and packages,
especially [`seriation`](https://github.com/mhahsler/seriation)
(Hahsler, Hornik, and Buchcta 2008) and
[`vegan`](https://CRAN.R-project.org/package=vegan) (Oksanen et al.
2024) that provide the means to seriate or ordinate matrices, especially
for frequency or count data. While binary (presence/absence) data are
often viewed as a reductive case of frequency data, they can also
present their own challenges. Moreover, not all incidence matrices (the
matrix of 0/1s that record the joint incidence or occurrence for a
row-column pairing) will necessarily be well seriated. The selection of
which row and column elements to inlcude in the input is accordingly an
intrinsic part of the task of seriation. In this respect, `lakhesis`
seeks to complement existing methods in `R`, focusing on binary data, by
providing an interactive, graphical means of selecting seriated
sequences. It uses correspondence analysis, a mainstay technique for
seriation, which is then fit to a reference curve that represents
“ideally” seriated data. Multiple seriations can be rerun on partial
subsets of the initial incidence matrix, which are then recompiled into
a single consensus seriation using an optimality criterion. The process
of harmonizing different partial strands of sequential elements via
iterative linear regressions is called a lakhesis technique, after the
fate from ancient Greek mythology figure who measured the strand of
one’s life. The package relies on `Rcpp` and `RcppArmadillo` for
Eddelbuettel and Balamuta (2018).

While command line functions can be run in `R`, the functionality of
`lakhesis` is primarily achieved via the Lakhesis Calculator, a
graphical platform in `shiny` (Chang et al. 2024) that enables
investigators to explore datasets for potential seriated sequences,
select them, and then harmonize them into a single consensus seriation.
Four panels are displayed in the calculator:

- **Seriation Explorer** (Top left) Displays the correspondence analysis
  of a dataset which has been fit to the curve an “ideal” seriation. Two
  plots are available, the biplot of the row and column CA scores as
  they have been fit to an ideally seriated curve, and a plot showing
  the orthogonal projection of the row and column scores as they have
  been fit to the reference curve. Selections can be made on either the
  biplot or the reference plot.
  - Map options. To the left, users can choose either a **symmetric** or
    an **asymmetric** plot.
- **Save Strand** (Top middle) Records the displayed plot as a partial
  seriation, or “strand” (i.e., partial with respect to the initial
  data).
  - Strands can be sequenced according to different projections:
    - **CA1** / **CA2** Projection of scores along the first or second
      principal CA axis.
    - **Procrustes1** / **Procrustes2** Projection of scores along the
      first or second axis after Procrustes fitting.
    - **Curve** Projection along the reference curve of an ideal
      seriation.
  - **Lakhesize** Produce a consensus seriation (must have saved at
    least two strands). Constructs a consensus seriation of the selected
    strands using an iterative process of linear regression of partial
    rankings in an agglomerative fashion. The matrix plot displays the
    incidence matrix of the resulting consensus seriation, with
    optimality criteria. The agreement of the seriation in each strand
    with that of the consensus seriation as well as its criterion is
    displayed in the Diganostics panel. The function `lakhesize()`
    performs this task.
  - **Run Deviance Test** Performs a goodness-of-fit test using
    deviance, treating the distribution of the row and column incidences
    with a quadratic-logistic model. The largest $p$ values of the row
    and column elements is contained in the Diganostics panel. The
    function `element_eval()` performs this task.
- **Consensus Seriation** (Top right) Displays the results of
  harmonizing selected partial seriations, which have been identified as
  “strands.” The process of deriving a consensus seriation entails a
  process of iterative regressions on partially seriated sequences,
  optimized using the concentration measure. The seriated incidence
  matrix is also displayed in this panel.
- **Diagnostics** (Bottom left) Critical coefficients to determine
  whether discordant strands should be removed and/or row or column
  elements should be suppressed from consideration.
  - **Agreement** expresses whether a strand agrees with consensus
    seriation.
  - **Criteria** expresses how well seriated the strand is. Options for
    optimality criteria are those which are used in Lakhesis,
    comprising:
    - **Squared correlation coefficient** (`cor_sq`)
    - **Weighted row-column concentration** (`conc_wrc`)
  - Tabs marked **Deviance** report on the goodness-of-fit of row and
    column elements in the consensus seriation using deviance with a
    quadratic-logistic model. Higher $p$ values will indicate poorer fit
    for a particular row or column element.
- **Modify** (Bottom right) Temporarily suppress row or column elements
  from correspondence analysis. Strands which have low agreement or high
  concentration may also be deleted in this panel.

The sidebar contains the following commands:

- **Choose CSV** Data must be without a header in a two-column “long”
  format of occurring pairs of row and column elements, where the first
  column contains a row element and the second column contains a column
  element of the incidence matrix.
- **Reinitialize** Resets the plots to their original, starting
  condition.
- **Replot with Selection** Upon the selection of row and column points
  from the Seriation Explorer panel, this command will perform and fit
  CA only on the selection. To return to the initial dataset, press the
  Reinitialize button. The function `ca_procrustes_ser()` performs this
  task.
- **Export Data** Will download results in a single `.rds` file, which
  is a `list` class object containing the following:
  - `consensus` The results of `lakhesize()`, a `lakhesis` class object
    containing row and column consensus seriations, coefficients of
    agreement and concentration, and the seriated incidence matrix.
  - `strands` The strands selected to produce `consensus`.

## Installation

To obtain the current development version of `lakhesis` from GitHub,
install from GitHub in the `R` command line with:

``` r
library(devtools)
install_github("scollinselliott/lakhesis", dependencies = TRUE, build_vignettes = TRUE) 
```

## Usage

To start the Lakhesis Calculator, execute the function `LC()`:

``` r
library(lakhesis)
LC()
```

In uploading a `csv` file for analysis inside the Lakhesis Calculator,
the incidence matrix should be in “long” format. That is, the file
should consist of just two columns without headers, in which each row
represents the incidence of a row-column pair. For example, an incidence
matrix of

$$\begin{array} \,  & C_1 & C_2 & C_3 \\\ R_1 & 1 & 0 & 0 \\\ R_2 &  0 & 1 & 1 \\\  R_3 & 0 & 0 & 1 \end{array}$$

will have a corresponding long format of

``` r
R1, C1
R2, C2 
R2, C3 
R3, C3
```

If characters are not displyaing properly in the plot, make sure to
check font encoding (UTF-8 is recommended).

Row and column elements must be unique (a row element cannot have the
same name as a column element).

The Lakhesis Calculator enables the temporary suppression of row or
column elements from the plots, with zero rows/columns automatically
removed. As such, unexpected results may be elicited if key elements are
suppressed. All elements can easily be re-added and the starting
incidence matrix re-initialized.

### Incidence Matrices

If data are already in incidence matrix format, the `im_long()` function
in `lakhesis` can be used to convert an incidence matrix to be exported
into the necessary long format, using the `write.table()` function to
export (see documentation on `im_long()`):

``` r
# x is a matrix of 0/1 values with unique row/column names
y <- im_long(x)
write.table(y, file = "im.csv", sep = ",")
```

The file `im.csv` can then be loaded into the Lakhesis Calculator.

## Consensus Seriations

Establishing a consensus seriation via a lakhesis technique can be done
in the calculator, but if one has seriations, whether derived by
Procustes-fit CA or by another method, one can perform a consensus
seriation in the console by creating a `strands` object and then
executing the `lakhesize()` function.

The console can also be used to perform consensus seriations. For
example, using the built-in selection of three strands in the data
object `qf_strands`, a consensus seriation is performed using the
`lakhesize()` function:

``` r
x <- lakhesize(qf_strands)
summary(x)
```

The vignette “A Guide to Lakhesis” contains more information on usage.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chang_shiny_2024" class="csl-entry">

Chang, W., J. Cheng, J. J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J.
Allen, J. McPherson, A. Dipert, and B. Borges. 2024. *Shiny: Web
Application Framework for R*. <https://shiny.posit.co>.

</div>

<div id="ref-eddelbuettel_extending_2018" class="csl-entry">

Eddelbuettel, D., and J. J. Balamuta. 2018. “Extending R with C++: A
Brief Introduction to Rcpp.” *The American Statistician* 72: 28–36.
<https://doi.org/10.1080/00031305.2017.1375990>.

</div>

<div id="ref-eddelbuettel_rcpparmadillo_2014" class="csl-entry">

Eddelbuettel, D., and C. Sanderson. 2014. “RcppArmadillo: Accelerating R
<span class="nocase">with</span>
<span class="nocase">high</span>-<span class="nocase">performance</span>
C++ <span class="nocase">linear</span>
<span class="nocase">algebra</span>.” *Computational Statistics and Data
Analysis* 71: 1054–63. <https://doi.org/10.1016/j.csda.2013.02.005>.

</div>

<div id="ref-hahsler_getting_2008" class="csl-entry">

Hahsler, M., K. Hornik, and C. Buchcta. 2008. “Getting Things in Order:
An Introduction to the R Package Seriation.” *Journal of Statistical
Software* 25: 1–34. <https://doi.org/10.18637/jss.v025.i03>.

</div>

<div id="ref-ihm_contribution_2005" class="csl-entry">

Ihm, P. 2005. “A Contribution to the History of Seriation in
Archaeology.” In *Classification – The Ubiquitous Challenge*, edited by
C. Weihs and W. Gaul, 307–16. Berlin: Springer.

</div>

<div id="ref-oksanen_vegan_2024" class="csl-entry">

Oksanen, J., G. L. Simpson, F. G Blanchet, R. Kindt, P. Legendre, P. R.
Minchin, R. B. O’Hara, et al. 2024. “Vegan: Community Ecology Package.”
<https://doi.org/10.32614/CRAN.package.vegan>.

</div>

<div id="ref-ter_braak_weighted_1986" class="csl-entry">

ter Braak, C. J. F., and C. W. N. Looman. 1986. “Weighted Averaging,
Logistic Regression and the Gaussian Response Model.” *Vegetatio* 65:
3–11. <https://doi.org/10.1007/BF00032121>.

</div>

</div>
