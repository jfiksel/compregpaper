source(file.path("..", "scripts", "compRegression.R"))
library(ggtern)
library(compositions)
library(robCompositions)
library(tidyverse)
source(file.path("..", "scripts", "multFracReg.R"))
data("WhiteCells", package = 'ggtern')
Image <- subset(WhiteCells, Experiment == "ImageAnalysis")
Microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")


ImageComp <- as.matrix(Image[,c("G", "L", "M")])
ImageComp  <- ImageComp  / rowSums(ImageComp )
MicroscopicComp <- as.matrix(Microscopic[,c("G", "L", "M")])
MicroscopicComp <- MicroscopicComp / rowSums(MicroscopicComp)

### Do CompReg Fit and get predicted lines
compregFit <- pairedCompReg(yout = MicroscopicComp, ypred = ImageComp)

### Plot coefficients
boot_out <- pairedCompReg(yout = MicroscopicComp, ypred = ImageComp, boot.ci = TRUE)
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



### Plot 
whitecells_coef_plot <-
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
ggsave("whitecellsCoefficients.pdf", whitecells_coef_plot,
       width = 8, height = 4)

### MFLR fit
XIlr <- ilr(ImageComp)
XIlrIntercept <- cbind(1, XIlr)
mflr_est <- fixedUpdateMFLR(MicroscopicComp, XIlrIntercept)

### ILR fit
ilrFit <- lm(ilr(MicroscopicComp) ~ ilr(ImageComp))




### Compare CompReg to other methods using LOO
ypred_compreg <- ypred_ilr <- ypred_mflr <- matrix(NA, nrow = nrow(MicroscopicComp), ncol = ncol(MicroscopicComp))
for(i in 1:nrow(MicroscopicComp)) {
    XIlr <- ilr(ImageComp)
    XIlrIntercept <- cbind(1, XIlr)
    ### CompRegFit
    compregFit <- pairedCompReg(yout = MicroscopicComp[-i,], ypred = ImageComp[-i,])
    ypred_compreg[i,] <- ImageComp[i,] %*% compregFit
    
    ### MFLR fit
    mflr_est <- fixedUpdateMFLR(MicroscopicComp[-i,], XIlrIntercept[-i,])
    ypred_mflr[i,] <- t(apply( XIlrIntercept[i,] %*% mflr_est,1,softmax))
    
    ### ILR fit
    ilrFit <- lm(ilr(MicroscopicComp[-i,]) ~ ilr(ImageComp[-i,]))
    ypred_ilr[i,] = ilrInv(XIlrIntercept[i,] %*% ilrFit$coefficients)
}


### Plot predicted fit
pred_plot_df <- data.frame(observed = rep(as.vector(MicroscopicComp), 3),
                           predicted = c(as.vector(ypred_compreg),
                                         as.vector(ypred_mflr),
                                         as.vector(ypred_ilr)),
                           method = rep(c('compreg', 'mflr', 'ilr'), each = nrow(MicroscopicComp) * ncol(MicroscopicComp)),
                           colname = rep(rep(c('Granulocytes', 'Lymphocytes', 'Monocytes'), each = nrow(MicroscopicComp)),3))
whitecells_pred_plot <-
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
ggsave("whitecellsPredictionsDirect.pdf", whitecells_pred_plot,
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

kldCompReg <- kldMat(MicroscopicComp, ypred_compreg)
kldILR <- kldMat(MicroscopicComp, ypred_ilr)
kldMFLR <- kldMat(MicroscopicComp, ypred_mflr)

mvarCompReg <- mvarCustom (MicroscopicComp, ypred_compreg)
mvarILR <- mvarCustom (MicroscopicComp, ypred_ilr)
mvarMFLR <- mvarCustom (MicroscopicComp, ypred_mflr)


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
    ll_null <- sum(yout * log(y_avg))
    
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
permutPval <- permutationCompReg(MicroscopicComp, ImageComp)

### ILR fit
ilrFit <- lm(ilr(MicroscopicComp) ~ ilr(ImageComp))
predFitILR = ilrInv(predict(ilrFit))


### Plot predicted fit
opar <- par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(3,3,0,0))
Y <- MicroscopicComp
for(i in 1:3){
    for(j in 1:3){
        plot(log(Y[,i]/Y[,j]),
             log(predFitCR[,i]/predFitCR[,j]),
             pch=ifelse(i!=j,1,""))
        if(i==j){
            text(x=0,y=0, labels=colnames(Y)[i],cex=1.5)
        }else{
            abline(a=0,b=1,col="gray",lwd=3)
        }
    }
}

colnames(Y) <- c("Granulocytes", "Lymphocytes", "Monocytes")
pdf("whitecellsPredictionsDirect.pdf", width = 8, height = 3)
opar <- par(mfrow=c(1,3), mar=c(2,2,1,1), oma=c(3,3,0,0))
for(i in 1:3){
    plot(Y[,i], predFitCR[,i], main = colnames(Y)[i], xlim=c(0,1), ylim=c(0,1))
    abline(a=0,b=1,col="gray",lwd=3)
}
mtext(text=c("Observed","Predicted"), side=c(1,2),at=0.5,line=1,outer=TRUE)
dev.off()