library(compositions)
library(robCompositions)
library(ggtern)
library(latex2exp)
library(cluster)
library(MASS)
library(tidyverse)
source("../scripts/compRegression.R")
source(file.path("..", "scripts", "multFracReg.R"))

data("educFM")
father <- as.matrix(educFM[,2:4])
father <- father / rowSums(father)
mother <- as.matrix(educFM[,5:7] )
mother <- mother/rowSums(mother)

### Exploratory plot--looks very linear
par(mfrow=c(3,3))
par(mar=rep(1, 4))
for(i in 1:3) {
    for(j in 1:3) {
        plot(mother[,j], father[,i], xlab="Mother", ylab = "Father")
    }
}

### Do CompReg Fit to get estimate of B
compregFit <- pairedCompReg(yout = father, ypred = mother)


### Plot coefficients
boot_out <- pairedCompReg(yout = father, ypred = mother, boot.ci = TRUE)
boot_out_mat <- boot_out$t
pp <- 0.95

ellipses_list <- lapply(1:3, function(i) {
    # find smallest ellipse containing the specified proportion of bootstrapped medians
    indices <- c(i, i + 3)
    fit <- cov.rob(boot_out_mat[,indices], quantile.used = ceiling(pp*nrow(boot_out_mat)), method = "mve")
    best_ellipse <- predict(ellipsoidhull(boot_out_mat[fit$best,indices] ))
    ellipse_df <- data.frame(y1 = best_ellipse[,1], y2 = best_ellipse[,2])
    ellipse_df[ellipse_df <0] <- 0
    ellipse_df[ellipse_df >1] <- 1
    ellipse_df$y3 <- pmax(0, 1 - ellipse_df$y1 - ellipse_df$y2)
    ellipse_df <- ellipse_df / rowSums(ellipse_df)
    return(ellipse_df)
})

est_df <- data.frame(y1 = compregFit[,1],
                     y2 = compregFit[,2],
                     y3 = compregFit[,3])

### Slight fudging so you can actually see the points
est_df[est_df < .001] <- .015

ellipses_list[[2]]$y1 <- ellipses_list[[2]]$y1 + .02

### Plot 
education_coef_plot <-
    ggtern(data = est_df , mapping = aes(x = y1, y = y2, z = y3)) +
    geom_polygon(data = ellipses_list[[1]], colour="blue", fill = NA) +
    geom_polygon(data = ellipses_list[[2]], colour="blue", fill = NA) +
    geom_polygon(data = ellipses_list[[3]], colour="blue", fill = NA)+
    scale_T_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    scale_L_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    scale_R_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    #geom_mean_ellipse(data = ellipses_list[[1]]) +
    #geom_mean_ellipse(data = ellipses_list[[2]]) +
   # geom_mean_ellipse(data = ellipses_list[[3]])+
    geom_text(aes(label = lab)) +
    xlab(TeX("$E\\[y_1\\]$")) +
    ylab(TeX("$E\\[y_2\\]$")) +
    zlab(TeX("$E\\[y_3\\]$")) 
ggsave("educationCoefficients.pdf", education_coef_plot,
       width = 8, height = 4)


### MFLR fit
XIlr <- ilr(mother)
XIlrIntercept <- cbind(1, XIlr)
mflr_est <- fixedUpdateMFLR(father, XIlrIntercept)

### ILR fit
ilrFit <- lm(ilr(father) ~ ilr(mother))

### Compare CompReg to other methods using LOO
ypred_compreg <- ypred_ilr <- ypred_mflr <- matrix(NA, nrow = nrow(father), ncol = ncol(father))
for(i in 1:nrow(father)) {
    XIlr <- ilr(mother)
    XIlrIntercept <- cbind(1, XIlr)
    ### CompRegFit
    compregFit <- pairedCompReg(yout = father[-i,], ypred = mother[-i,])
    ypred_compreg[i,] <- mother[i,] %*% compregFit
    
    ### MFLR fit
    mflr_est <- fixedUpdateMFLR(father[-i,], XIlrIntercept[-i,])
    ypred_mflr[i,] <- t(apply( XIlrIntercept[i,] %*% mflr_est,1,softmax))
    
    ### ILR fit
    ilrFit <- lm(ilr(father[-i,]) ~ ilr(mother[-i,]))
    ypred_ilr[i,] = ilrInv(XIlrIntercept[i,] %*% ilrFit$coefficients)
}


### Plot predicted fit
pred_plot_df <- data.frame(observed = rep(as.vector(father), 3),
                           predicted = c(as.vector(ypred_compreg),
                                         as.vector(ypred_mflr),
                                         as.vector(ypred_ilr)),
                           method = rep(c('compreg', 'mflr', 'ilr'), each = nrow(father) * ncol(father)),
                           colname = rep(rep(c('Low', 'Medium', 'High'), each = nrow(father)),3))
pred_plot_df$colname <- factor(pred_plot_df$colname, levels = c('Low', 'Medium', 'High' ))

edfm_pred_plot <-
    pred_plot_df %>%
    filter(method == 'compreg') %>%
    ggplot(aes(x = observed, y = predicted)) +
    geom_point(size = 2, shape = 1) +
    geom_abline(alpha = .4) +
    facet_wrap(~colname) +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw() +
    xlab("Observed") +
    ylab("Predicted")
ggsave("educationPredictionsDirect.pdf", edfm_pred_plot,
        width = 8, height = 4)



### KLD
kld <- function(true, est) {
    d <- 0
    for(i in 1:length(true)) {
        if(true[i] == 0) {
            d <- d 
        } else {
            d <- d + true[i] * log(true[i] / est[i])
        }
    }
    return(d)
}

kldMat <- function(trueMat, estMat) {
    d <- 0
    N <- nrow(trueMat)
    for(i in 1:N) {
        d <- d + kld(trueMat[i,], estMat[i,])
    }
    return(d / N)
}

mvarCustom <- function(truth, pred) {
    1 - mvar(acomp(truth) - acomp(pred)) / mvar(acomp(truth))
}

kldCompReg <- kldMat(father, ypred_compreg)
kldILR <- kldMat(father, ypred_ilr)
kldMFLR <- kldMat(father, ypred_mflr)

mvarCompReg <- mvarCustom (father, ypred_compreg)
mvarILR <- mvarCustom (father, ypred_ilr)
mvarMFLR <- mvarCustom (father, ypred_mflr)

### 



### Permutation test
logLikComp <- function(y, x, M) {
    mu <- x %*% M
    return(sum(y*log(mu)))
}


permutationCompReg <- function(yout, ypred, nperms = 1000, init.seed = NULL) {
    if(is.null(init.seed)){
        set.seed(123)
    } else {
        set.seed(init.seed)
    }
    ### Get observed M
    m_obs <- pairedCompReg(yout, ypred)
    ### Get observed log-likelihood
    ll_full <- logLikComp(yout, ypred, m_obs)
    ### Get observed log-likelihood under null
    y_avg <- colMeans(yout)
    ll_null <- sum(sweep(yout, MARGIN = 2, log(y_avg), `*`))
    
    llr_obs <- ll_full - ll_null
    ### Get do permutation test
    permut_stats <- rep(NA, nperms)
    for(i in 1:nperms) {
        perm_index <- sample(1:nrow(ypred))
        ypred_perm <- ypred[perm_index,]
        m_perm <- pairedCompReg(yout, ypred_perm)
        ll_perm <- logLikComp(yout, ypred_perm, m_perm)
        permut_stats[i] <- ll_perm - ll_null
    }
    pval <- mean(permut_stats >= llr_obs)
    return(pval)
}
permutPval <- permutationCompReg(father, mother)

