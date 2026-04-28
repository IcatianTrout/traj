# traj 3.0.2

-   Fixed a small bug that occurred in `trajMeasures` when ID=FALSE and some trajectories have less than 3 observations.
-   Changed the `trajdata` data set by making trajectories 1 and 130 into outliers.
-   Changed the I


# traj 3.0.1

-  In the `plot` function, setting the argument `sample.size` to `NULL` plots all the trajectories on the same graph (in a random order).
-   In the `trajClusters` function, setting the argument `subset.n` to a numerical integer while `nclusters` is set to `NULL` makes it so the optimal number of clusters is determined using a random sample of the data of size `subset.n`. Assuming `subset.n` is large enough that the random sample is representative of the data, this would speed up the process of identifying the optimal number of clusters.
-   In `scatterplots`, the legend now appear outside the scatter plots, for improved visibility.
-   In the `scatterplots` function, we added an argument `which.scatter` allowing to plot only a subset of all the available scatter plots.
-   In the `scatterplots` function, we added an argument `N` allowing to plot a random sample of size `N`, while preserving the groups' relative sizes. Assuming `N` is large enough that the sample is representative of the data, this would speed up the plotting process.


# traj 3.0.0

-   Substantial modifications to the measures.
-   The clustering step now relies on a version of the Spectral Clustering algorithm.
-   The main function are `trajMeasures` (computes the measures) and `trajClusters` (finds the clusters), with `trajReduce` (finds a representative subset of measures) being accessory.
-   The plotting functions that `trajClusters` can be passed into are now plot, `scatterplots` and `CVIplot`.
-   The `trajClusters` function responsible of finding the clusters has a logical argument that allows to choose between soft and hard clustering.


# traj 2.2.1

-   Added a new measure, "m5: slope of linear model".
-   Spit the "plot" function into three: `plot`, `plot.scatter` and `plot.crit`.
-   Improved the presentation of the scatter plots.
-   Made `Step1Measures` more lenient with how the input data is formatted.

# traj 2.2.0

-   Added k-medoids as the default clustering algorithm.
-   Added the Calinski-Harabasz index as the default criterion for determining the optimal number of clusters.

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
