# traj 2.1.0

-   Makes substantial changes to the list of measures.
-   In `trajdata`, the group of size 30 is made up of quadratic (instead of linear) curves.
-   Introduces the vignette "Using the traj package".

# traj 2.0.1

-   Fixes minor bugs and improves presentation.

# traj 2.0.0

-   Changes the way the measures are computed in order to be better suited to missing values and unequally spaced observation times.
-   Allows for better control over the choice of a "midpoint".
-   Implements a capping procedure to automatically handle outliers. Measures of the form 0/0 are set to 1.
-   Fixes an important bug in the implementation of step 2.
-   The criterion of Tibshirani et al. based on the GAP statistic is the new default for choosing the number of clusters.
-   Allows more control over the parameters of `kmeans`.
-   The outputs of `summary`, `print` and `plot` are more detailed.
