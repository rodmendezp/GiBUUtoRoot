# Using Inverse Probability Weighting with SVMs to Address Confounding
### by Kristin A. Linn
### June 25, 2015

Here we provide an example of how to implement inverse probability weighting with SVMs to address confounding.  The basic setup is that we have feature, class label pairs of the form (
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_03.png" alt="Equation Fail"height="20">
) for each subject, i=1,...,n, where 
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_04.png" alt="Equation Fail"height="20">
 and 
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_05.png" alt="Equation Fail"height="20">
 for all 
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_06.png" alt="Equation Fail"height="20">
. We wish to train a SVM to predict y given x. As an example, y might be an indicator of disease/control group and x might be a vectorized image containing voxel values or volumes of regions across the brain. However, the additional feature vector 
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_07.png" alt="Equation Fail"height="20">
 observed for all subjects confounds the relationship between x and y. For example, a might contain covariates such as age and sex.  In the presence of confounding by a, inverse probability weighting is used to recover an estimate of the target classifier, which is the SVM classifier that would have been estimated had there been no confounding by a.



```r
set.seed (1)
```

We use the package 'rPython' to access libSVM (https://www.csie.ntu.edu.tw/~cjlin/libsvm/) through scikit learn (http://scikit-learn.org/stable/). The file fit_svm.py contains a python function that implements a linear kernel SVM with subject-level weights and a grid search to tune the cost parameter, C.



```r
print("I am R Code")
```

```
## [1] "I am R Code"
```

## Generate data for example

We generate data such that the confounders, a1 and a2, affect both the features, x1 and x2, as well as the class labels, d.


N.B. If some of the estimated weights are extremely large, one may consider truncating the predicted probabilities (e.g., at the 5th percentile) or using stabilized weights. Define 
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_08.png" alt="Equation Fail"height="20">
 and 
<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_09.png" alt="Equation Fail"height="20">
 as well as corresponding estimates \hat{S}_i = \hat{pr}(d_i=1 | a1_i, a2_i) and \hat{M}_i = \hat{pr}(d_i=1), Then, stabilized weights and their corresponding estimates are defined, respectively, as:

<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_01.png" alt="Equation Fail"height="20">

and fun

<img src="https://github.com/rodmendezp/GiBUUtoRoot/blob/master/Particles/eq_no_02.png" alt="Equation Fail"height="20">

## Train the inverse probability weighted SVM (IPW-SVM)


