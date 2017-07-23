# JDINAC: joint density-based non-parametric differential interaction network analysis and classification using high-dimensional sparse omics data
<br >

**Cite:** Ji, J., He, D., Feng, Y., He, Y., Xue, F., & Xie, L. (2017). JDINAC: joint density-based non-parametric differential interaction network analysis and classification using high-dimensional sparse omics data. _Bioinformatics_. [doi: 10.1093/bioinformatics/btx360](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btx360)

<br >

## Usage

R code avaiable from https://github.com/jijiadong/JDINAC/blob/master/jdinac.r
```
jdinac(EDGE, classLabel, DataFit, DataPre, zFit=NULL, zPre=NULL, nsplit=10, nfolds=5)
```
### Input:

  - DataFit, DataPre: data matrices containing one sample per row, one variable per column. DataFit: training data, DataPre: prediction data. 
  - zFit, zPre: covariate in training data and prediction data respectively.
  - classLabel: must be 0 or 1, e.g. 1 for cases and 0 for controls.
  - nsplit: randomly split the dataset *2\*nsplit* times.
  - nfolds: number of folds for Cross-validation.
  - EDGE: array indices. which edge will be tested.

### Output

  - yPre: predicted value of response.
  - Eset: differential edge set; the first two columns are the array indices; the edge between *row*-th variable and *col*-th variable. The 3rd column (*numb*) is the differential dependency weight; not normalized.

### Example
Run the R code *jdinac.r* first. 

```{r}
library(mvtnorm)
Sigma0 <- 0.5* (-1)^abs(outer(1:10,1:10,"-"))
Sigma1 <- array(0.5,dim=c(10,10))
adjM <- array(1,dim=c(10,10))
adjM[lower.tri(adjM,diag=T)] <- 0
adjM[1,5] <- adjM[3,6] <- adjM[4,8] <- 0
adjM
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
##  [1,]    0    1    1    1    0    1    1    1    1     1
##  [2,]    0    0    1    1    1    1    1    1    1     1
##  [3,]    0    0    0    1    1    0    1    1    1     1
##  [4,]    0    0    0    0    1    1    1    0    1     1
##  [5,]    0    0    0    0    0    1    1    1    1     1
##  [6,]    0    0    0    0    0    0    1    1    1     1
##  [7,]    0    0    0    0    0    0    0    1    1     1
##  [8,]    0    0    0    0    0    0    0    0    1     1
##  [9,]    0    0    0    0    0    0    0    0    0     1
## [10,]    0    0    0    0    0    0    0    0    0     0

EDGE <- which(adjM!=0, arr.ind=T)
head(EDGE)
##      row col
## [1,]   1   2
## [2,]   1   3
## [3,]   2   3
## [4,]   1   4
## [5,]   2   4
## [6,]   3   4

set.seed(2017)
size0 <- size1 <- 50
class0 <- rmvnorm( n = size0, sigma = Sigma0,method = "svd")
class1 <- rmvnorm( n = size1, sigma = Sigma1,method = "svd" ) 
dataset <- rbind(class0,class1)
classLabel <- c(rep(0,size0),rep(1,size1))

difnet <- jdinac(EDGE=EDGE,classLabel=classLabel,DataFit=dataset,DataPre=dataset,nsplit=10,nfolds=5)
head(difnet$yPre)
## [1] 3.893804e-01 3.459244e-08 9.624137e-07 8.148533e-09 4.774174e-07
## [6] 3.384799e-05

head(difnet$Eset)
##    row col numb
## 1    1   2   20
## 3    2   3   16
## 5    2   4   14
## 6    3   4   13
## 11   2   6   13
## 21   2   8   13

```

