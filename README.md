# MSA_FA_LDA

**Project Description:**
Performing Factor Analyses technique for dimensionality reduction and Linear Discriminant Analyses on different dataset for a university group project (working in R)


**Assignment 1 (FA):**
The data set pulp_paper contains measurements of properties of pulp bers and the paper made from them. There are n = 62 observations on 4 paper properties: (BL) breaking length, (EM) elastic modulus, (SF) stress at failure, (BS) burst strength; and 4 pulp ber characteristics: (AFL) arithmetic ber length, (LFF) long ber fraction, (FFF) ne ber fraction, (ZST) zero span tensile.

1) Obtain the maximum likelihood solution for $m = 2$ and $m = 3$ common factors on the standardize observations and compute the proportion of total sample variance due to each factor. List the estimated communalities, specic variances, and the residual matrix. Compare the results. Which choice of m do you prefer? Why?
2) Give an interpretation to the common factors in the $m = 2$ solution.
3) Make a scatterplot of the factor scores for $m = 2$ obtained by the regression4method. Is their correlation equal to zero? Should we expect so? Comment.
4) Suppose we have a new observation (15.5; 5.5; 2; -0.55; 0.6; 65; -5; 1.2). Calculate the corresponding $m = 2$ factor scores and add this bivariate point toù the plot in 3). How is it placed compared to the rest of the $n = 62$ points?
Could you tell without computing the factor scores? Comment.


**Assignment 2 (LDA):**
The dataset glass contains data on n = 214 single glass fragments. Each case has a measured refractive index (RI) and composition (weight percent of oxides of Na, Mg, Al, Si, K, Ca, Ba and Fe). The composition sums to around $100%$; what is not anything else is sand. The fragments are classied as six types (variable type). The classes are window  float glass (WinF), window non  float glass (WinNF), vehicle window glass (Veh), containers (Con), tableware (Tabl) and vehicle headlamps (Head).

1) Use linear discriminant analysis to predict the glass type. Look at the first two discriminant directions: what are the most important variables in separating the classes? Comment.
2) Compute the training error. Are there any groups less homogeneous than the others? Comment.
3) Implement a 10-fold cross validation using the partition of the observations provided by the variable groupCV to estimate the error rate. Comment.
4) Use the first two discriminant variables for a two-dimensional representation of the data together with centroids by using color-coding for the 6 classes of the class variable type. Comment in view of the answer to point 2).
5) Compute the training error and the 10-fold cross validation error for each reduced-rank LDA classier. Plot both error curves against the number of discriminant directions, add full-rank LDA errors found in points 2) and 3). What classifier do you prefer? Comment.
