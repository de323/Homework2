# Homework2

    ridge_regression
    @title Compute Ridge Regression Coefficient Vector
    @param X A numeric data matrix
    @param y Response vector.
    @param lambda_vals A sequence of penalty terms

    @examples X <- as.matrix(iris[,-c(3,5)])
             Y <- iris[[3]]
             ridge_regression(X=X, y=Y, lambda_vals = 10
    @references
    Arnold, T. (2019). Statsmaths/casl [R]. Retrieved from https://github.com/statsmaths/casl (Original work published 2018)
    Arnold, T., Kane, M., & Bryan, L. (2019). A Computational Approach to Statistical Learning.
    How and when: Ridge regression with glmnet • blogR. (n.d.). Retrieved October 23, 2019, from BlogR on Svbtle website: 
    https://drsimonj.svbtle.com/ridge-regression-with-glmnet
    Introduction to roxygen2. (n.d.). Retrieved October 23, 2019, from https://cran.r-
    project.org/web/packages/roxygen2/vignettes/roxygen2.html
    R—Why is glmnet ridge regression giving me a different answer than manual calculation? - Cross Validated. (n.d.).     
    Retrieved October 23, 2019, from https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-
    me-a-different-  answer-than-manual-calculat
    (N.d.). Retrieved October 23, 2019, from https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
    Collaborations with: Moid Ali and Tyler Harvey
    Additional support from SAS/R tutor through the Yale School of Public Health.

    ridge_reg <- function(X, y, lambda_vals){
    svd_obj <- svd(X)
    U <- svd_obj$u
    V <- svd_obj$v
    svals <- svd_obj$d
    k <- length(lambda_vals)

    ridge_beta <- matrix(NA_real_, nrow = k, ncol = ncol(X))
    for (j in seq_len(k))
    {
      D <- diag(svals / (svals^2 + lambda_vals[j]))
      ridge_beta[j,] <- V %*% D %*% t(U) %*% y
    }

    ridge_beta

    return(list(coefficients=ridge_beta))
    }
    
    Cross Validation Lambda
    @param y A response variable.
    @param x A numeric data matrix
    @param nfolds Specifies the number of folds in the model

    @return The optimal lambda associatied with the least mean squares prediction error.

    @examples lambdaval=10
           mod <- lm_ridge_2(X, y, nfolds=10, lambda=lambdaval)
           mod$coefficients
    @references
    Arnold, T. (2019). Statsmaths/casl [R]. Retrieved from https://github.com/statsmaths/casl (Original work published 2018)
    Arnold, T., Kane, M., & Bryan, L. (2019). A Computational Approach to Statistical Learning.
    How and when: Ridge regression with glmnet • blogR. (n.d.). Retrieved October 23, 2019, from BlogR on Svbtle website: 
    https://drsimonj.svbtle.com/ridge-regression-with-glmnet
    Introduction to roxygen2. (n.d.). Retrieved October 23, 2019, from https://cran.r-
    project.org/web/packages/roxygen2/vignettes/roxygen2.html
    R—Why is glmnet ridge regression giving me a different answer than manual calculation? - Cross Validated. (n.d.). 
    Retrieved October 23, 2019, from https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-
    me-a-different-answer-than-manual-calculat
    (N.d.). Retrieved October 23, 2019, from https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
    Collaborations with: Moid Ali and Tyler Harvey
    Additional support from SAS/R tutor through the Yale School of Public Health.

    lm_ridge_2<- function (y,x, nfolds){
    n <- length(y)
    k <- length(lambda)
    di <- dim(y)[2]
    p <- dim(x)[2]
    msp <- matrix(nrow = nfolds, ncol = k)

    for (j in 1:nfolds) {
    y_test <- y[ nfolds[[ j ]], ]
    y_train <- y[ -nfolds[[ j ]], ]

    my <- lapply(y_train, mean, na.RM=TRUE)
    yy <- as.matrix( y_train + my)
    x_train <- x[ -nfolds[[ j ]] ]
    x_test <- x[ nfolds[[ j ]]]

    sa <- svd(xtrain)
    d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)
    d2 <- d^2    ;    A <- d * tu %*% yy

    for (i in 1:k) {
    beta <- crossprod(v / (d2 + lambda[i] ), A)
    est <- xtest %*% beta
    mode(est)="numeric"
    est_new <- data.frame(est) /my
    msp <- lapply((ytest - est_new)^2, mean, na.rm=TRUE)

    return(list(coefficients=lambda_vals[which.min(msp)]))
    }}}
    
    Tests for Ridge Regression & Cross Validation
    library(glmnet)
    library(testthat)
    usethis::use_testthat()

    context("Test the output of homework 2.")

    test_that("Your ridge regression() function works with iris data.", {

    X <- as.matrix(iris[,-c(3,5)])
    Y <- iris[[3]]
    n <- nrow(X)

    sd_y <- sqrt(var(Y)*(n-1)/n)

    fit_glmnet <- glmnet(X,Y, alpha=0, standardize= FALSE, intercept= FALSE, thresh= 1e-20)

    beta_new <- as.vector(coef(fit_glmnet, s=sd_y*10/n, exact= TRUE, x=X, y=Y))[-1]

    fit_ridge_reg <- ridge_reg(X=X, y=Y, lambda_vals= 10)

    expect_equivalent(beta_new, fit_ridge_reg$coefficients, tolerance = 1e-5)
    })

    test_that("Your ridge_regression() passes a test with a randomly generated dataframe.", {

    set.seed(123)
    n <- 1000
    p <- 100

    X <- matrix(rnorm(n*p, 0, 1), n,p)
    beta <- rnorm(p, 0, 1)
    Y <- X %*% beta + rnorm(n, 0, 0.5)

    sd_y <- sqrt(var(Y)*(n-1)/n)

    fit_glmnet <- glmnet(X, Y, alpha =0, standardize= FALSE, intercept= FALSE, thresh= 1e-20)

    beta_new <- as.vector(coef(fit_glmnet, s=sd_y*10/n, exact= TRUE, x=X, y=Y))[-1]

    fit_ridge_reg <- ridge_reg(X=X, y=Y, lambda_vals= 10)

    expect_equivalent(beta_new, fit_ridge_reg$coefficients, tolerance= 1e-5)
    })

    test_that("Your ridge_regression() passes cross validation with iris data.", {
    y <- as.numeric(iris[,5])
    X <- iris[y!=1, 1:4]
    y <- y[y!=1]-2

    n_sample = NROW(X)

    w = .6

    X_train = X[0:(w * n_sample),]
    y_train = y[0:(w * n_sample)]

    X_test = X[((w * n_sample)+1):n_sample,]
    y_test = y[((w * n_sample)+1):n_sample]

    set.seed(0)
    model_lambda <- cv.glmnet(as.matrix(X_train), as.factor(y_train), nfolds = 10,
                alpha=0, family="binomial", type.measure="class")

    best_lambda  <- model_lambda$lambda.1se
    best_lambda

    lambdaval=10
    mod <- lm_ridge_2(X, y, nfolds=10, lambda=lambdaval)
    mod$coefficients

    expect_lt(mo$coefficients-best_lambda, .001)
    })
    
    ---
    title: "Ridge-Regression"
    author: "Diana Estefania Estrada Alamo"
    date: "10/23/2019"
    output:
      html_document: default
      pdf_document: default
    ---

    ### 1.) **CASL 2.11 Exercises problem number 5. Consider the simple regression model with only a scalar x and intercept:**
    \[
    y = \beta_o \space + \space \beta_1 \times x
    \]

    ### **Using the explicit formula for the inverse of a 2-by-2 matrix, write down the least squares estimators for β_0^ and     
    β_1^**

    \[
    \hat{\beta} = ({X^T} X)^{-1}{X^T}Y, \space \space \space X \in \mathbb{R}^{n\times 2}
    \]

    - #### a. Applying the general formula 

    \[
    where \space A^{-1}, \space \space \space A \in \mathbb{R}^{2\times 2}
    \]

    \[
    A = 
    \left[\begin{array}
    {rrr}
    a & b\\
    c & d\\
    \end{array}\right], \space \space \space
    A^{-1} = \frac{1}{ad-bc}
    \left[\begin{array}
    {rrr}
    d & -b\\
    -c & a\\
    \end{array}\right]
    \]

    - #### b. Answer

    \[
    where \space ({X^T} X)^{-1} {X^T}Y, 
    \space \space \space \mathbb{R}^{2 \times 1}, 
    \space \space \space \beta = 
    \left[\begin{array}
    {rrr}
    \beta_0\\
    \beta_1\\
    \end{array}\right]
    \]

    \[
    {X^TX} = 
    \left[\begin{array}
    {rrr}
    \space 1 ... 1 \space \space\\
    X_1...X_n \\
    \end{array}\right] \times
    \left[\begin{array}
    {rrr}
    1 \space \space \space X_1\\
    \space : \space \space : \space\\
    1 \space \space \space X_n
    \end{array}\right] = 
    \left[\begin{array}
    {rrr}
    n &\sum{X_i} \space \space\\
    \sum{X_i} & \sum{X_i}^2 \\
    \end{array}\right]
    \]

    \[
    ({X^T}X)^{-1} = \frac{1}{n \sum {{X_j}^2}-(\sum {X_i})^2} \space
    \left[\begin{array}
    {rrr}
    \sum{X_i}^2 &-\sum{X_i}\\
    -\sum{X_i} & n\\
    \end{array}\right]
    \]

    \[
    {X^T}Y = 
    \left[\begin{array}
    {rrr}
    \space 1 ... 1 \space \space\\
    X_1...X_n \\
    \end{array}\right]
    \left[\begin{array}
    {rrr}
    Y_1\\
    \space \space : \space\\
    Y_n
    \end{array} \right] =
    \left[\begin{array}
    {rrr}
    \space \sum{Y_i} \space \space \\
    \sum{X_iY_i}\\
    \end{array}\right]
    \]

    \[
    2X^{T} X \hat{\beta} = 2X^{T} Y
    \]

    \[
    ({X^T} X)^{-1} {X^T} X \hat{\beta} = (X^{T} X)^{-1}X^{T} Y
    \]
    
    \[
    \hat{\beta} =  ({X^T} X)^{-1} {X^T}Y
    \]

    \[
    \hat{\beta} = \frac{1}{n \sum {{X_j}^2}-(\sum {X_i})^2} \space 
    \left[\begin{array}
    {rrr}
    \sum{X_i}^2 & -\sum{X_i}\\
    -\sum{X_i} & n\\
    \end{array}\right] \times
    \left[\begin{array}
    {rrr}
    \space \sum{Y_i} \space \space\\
    \sum{X_iY_i}\\
    \end{array}\right] 
    \]

    \[
    = 
    \frac{1}{n \sum {{X_j}^2}-(\sum {X_i})^2} \space 
    \left[\begin{array}
    {rrr}
    \sum{X_i}^2 \sum{Y_i} \space - \space \sum{X_i} \sum{X_iY_i} \\
    -\sum{X_i} \sum{Y_i} \space + \space n\sum{X_iY_i} \space \space \space \space\\
    \end{array}\right]
    \]

    \[If \space we \space suppose \space that \space X^{T}X \space is \space invertible, \space then \space the \space    
    estimated \space coefficient \space \hat{\beta} \space is \space a \space fixed \space \\ 
    linear \space combination \space of \space Y \space derived \space by \space X^{T}X^{-1}X^{T} \space as \space seen \space 
    above. \space \\
    We \space are \space trying \space to \space minimize \space the\space sum \space of \space square \space residuals \space    
    over \space our \space choice \space parameter \space vector \space \hat{\beta}, \\ 
    differentiating \space with \space respect \space to \space a \space vector \space with \space p \space independent \     
    paramenters.
    \]

    ###  4.) **Section 2.8 of CASL shows that as the numerical stability decreases, statistical errors increase. Reproduce the    
    results and then show that using ridge regression can increase numerical stability and decrease statistical error.**

    \[
    Av \space + \space e = z
    \]

    \[
    \lVert Av - A \hat{v} \lVert \space \le \space \epsilon
    \]

    \[
    \lVert Av - A \hat{v} \lVert \space = \space \lVert (Av - z) \space + \space (z \space - \space A \hat{v}) \lVert
    \]

    \[
    \le \space \lVert Av - z \lVert \space + \space \lVert z \space - \space A \hat{v} \lVert
    \]

    \[
    = \lVert e \lVert \space + \space \space \lVert Av - z \lVert
    \]

    \[
    \le \lVert e \lVert \space + \space \epsilon
    \]

    \[
    \frac {\lVert \hat{v} \space - \space v \lVert}{\lVert v \lVert \space} \le \frac {\lVert A(v -\hat{v}) \lVert}{\lVert Av       
    \lVert \space} 
    \]

    \[
    \frac {\lVert \hat{v} \space - \space v \lVert}{\lVert v \lVert \space} \le \frac {1}{\kappa_2(A) \times \lVert Av \lVert} 
    \times (\lVert e \lVert \space + \space \epsilon) 
    \]

    ```{r }
    devtools::install_github("statsmaths/casl")
    library(casl)

    n <- 1000; p <- 25
    beta <- c(1, rep(0, p - 1))
    X <- matrix(rnorm(n * p), ncol = p)

    svals <- svd(X)$d
    max(svals) / min(svals)

    N <- 1e4; l2_errors <- rep(0, N)
    for (k in 1:N) {
    y <- X %*% beta + rnorm(n)
    betahat <- casl_ols_svd(X, y)
    l2_errors[k] <- sqrt(sum((betahat - beta)^2))
    }
    mean(l2_errors)

    alpha <- 0.001
    X[,1] <- X[,1] * alpha + X[,2] * (1 - alpha)
    svals <- svd(X)$d
    max(svals) / min(svals)

    N <- 1e4; l2_errors <- rep(0, N)
    for (k in 1:N) {
    y <- X %*% beta + rnorm(n)
    betahat <- solve(crossprod(X), crossprod(X, y))
    l2_errors[k] <- sqrt(sum((betahat - beta)^2))
    }
    mean(l2_errors)
    ```
    
    \[
    As \space shown \space with \space the \space above \space evidence, \space the \space condition \space number \space is    
    \space an \space important \space player \space in \space numerical \space stability. \space \\
    Specifically, \space numerical \space instability \space highlights \space noise \space from \space arithemetic \space 
    floating \space point \space \\ 
    and \space statistical \space error \space from \space noise \space in \space data. \space These \space have \space an  
    \space inverse \space relationship \space with \space one \space another.
    \]

    ### 5.) **Consider the LASSO penalty. Show that if | X_j^T Y | ≤ nλ, then β_lasso must be zero.**

    \[
    L=\frac{1}{2n} \space \lVert \space {Y} - {X\beta} \space \lVert^2_2  \space + \space {\lambda} \space \lVert \space  
    {\beta} \space \lVert_1, \space Y\in \mathbb{R}^{n\times 1}, \space X\in \mathbb{R}^{n\times p}, \space \beta \in 
    \mathbb{R}^{p \times 1} 
    \]

    \[
    |\space {X_j}^T{Y} \space| < n \lambda 
    \]

    - #### a. Assumptions 

    \[
    ({X}^T{X}) = I, \space X_j...X_p \space are \space independent 
    \]

    \[
    {X_j}^TY-n \lambda < 0, \space \hat{\beta}_{lasso_j} = 0 
    \]

    \[
    \frac{1}{2n} \space \lVert \space {X} - {Y\beta} \space \lVert^2_2  \space + \space {\lambda} {\beta}, \space \beta > 0
    \]

    - #### b. Answer

    \[
    \frac{\delta L}{\delta \beta} = \frac{2}{2n} (-Y^T)(X - Y \beta) \space + \space \lambda = 0 
    \]

    \[
    = (-Y^T)(X - Y \beta) \space + \space \lambda = 0 
    \]

    \[
    = -Y^TX \space + Y^T Y \beta \space + n \lambda = 0 
    \]

    \[
    = Y^T Y \beta = Y^T X - n \lambda 
    \]

    \[
    \beta = (Y^T Y)^{-1}[Y^T X - n \lambda] 
    \]

    \[
    Y^T X - n \lambda = 
    \left[\begin{array}
    {rrr}
    ...Y_1... \\
    ...Y_0... \\
    ...Y_p...
    \end{array}\right] 
    X-n 
    \left[\begin{array}
    {rrr}
    \lambda_1\\
    \space \space : \space\\
    \lambda_p
    \end{array}\right] =
    \left[\begin{array}
    {rrr}
    {Y_i}^T X - n \lambda\\
    \space : \space \space \space \space \space \space \space \\
    {Y_p}^T X - n \lambda
    \end{array}\right] 
    \]

    ## **Citations** 

    #### An Example R Markdown. (n.d.-a). Retrieved October 23, 2019, from    
    http://www.math.mcgill.ca/yyang/regression/RMarkdown/example.html
    #### An Example R Markdown. (n.d.-b). Retrieved October 23, 2019, from 
    http://www.statpower.net/Content/312/Lecture%20Slides/SampleMarkdown.html
    #### How It Works. (n.d.). Retrieved October 23, 2019, from https://rmarkdown.rstudio.com/lesson-2.html
    #### Markdown Basics. (n.d.). Retrieved October 23, 2019, from https://rmarkdown.rstudio.com/authoring_basics.html
    #### math mode—How to write Euclidean distance. (n.d.). Retrieved October 23, 2019, from TeX - LaTeX Stack Exchange 
    website: https://tex.stackexchange.com/questions/70020/how-to-write-euclidean-distance
    #### Mathematics in R Markdown. (n.d.). Retrieved October 23, 2019, from 
    https://www.calvin.edu/~rpruim/courses/s341/S17/from-class/MathinRmd.html
    #### Pimp my RMD: a few tips for R Markdown. (n.d.). Retrieved October 23, 2019, from https://holtzy.github.io/Pimp-my-
    rmd/
    #### Rmarkdown-cheatsheet-2.0. (n.d.). 2.
    #### Rmarkdown-reference. (2014). 5.
    #### Using R Markdown for Class Reports. (n.d.). Retrieved October 23, 2019, from 
    https://www.stat.cmu.edu/~cshalizi/rmarkdown/#mark-up-markdown
    #### Arnold, T., Kane, M., & Bryan, L. (2019). A Computational Approach to Statistical Learning.
    #### YouTube. (n.d.). Retrieved October 23, 2019, from https://www.youtube.com/watch?v=C-uW45FSsNQ
