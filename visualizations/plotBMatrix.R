library(ggtern)
library(latex2exp)
library(ggpubr)
library(gtools)
library(cluster)
library(MASS)
source(file.path("..", "scripts", "compRegression.R"))

B1_df <- data.frame(y1 = c(.9, .05, .05),
                    y2 = c(.05, .9, .05),
                    y3 = c(.05, .05, .9),
                    x = 1:3)
B2_df <- data.frame(y1 = c(.4, .3, .3),
                    y2 = c(.3, .4, .3),
                    y3 = c(.3, .3, .4))
lab = 1:3
B1_plot <-
    ggtern(data = B1_df, mapping = aes(x = y1, y = y2, z = y3)) +
    #geom_point()  +
    geom_text(aes(label = lab)) +
    xlab(TeX("$E\\[y_1\\]$")) +
    ylab(TeX("$E\\[y_2\\]$")) +
    zlab(TeX("$E\\[y_3\\]$")) +
    ggtitle(TeX("$\\mathbf{B}^{1}")) +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_T_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    scale_L_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    scale_R_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) 
     

B2_plot <-
    ggtern(data = B2_df, mapping = aes(x = y1, y = y2, z = y3)) +
    #geom_point()  +
    geom_text(aes(label = lab)) +
    xlab(TeX("$E\\[y_1\\]$")) +
    ylab(TeX("$E\\[y_2\\]$")) +
    zlab(TeX("$E\\[y_3\\]$")) +
    ggtitle(TeX("$\\mathbf{B}^{2}")) +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_T_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    scale_L_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) +
    scale_R_continuous(labels = seq(.2, 1, by =.2), breaks = seq(.2, 1, by = .2)) 

B_viz_plot <- ggtern::grid.arrange(B1_plot, B2_plot, nrow =1)
ggsave("B_viz.pdf", width = 6, height = 3, plot = B_viz_plot)


#### Bootstrap
dir_data_gen <- function(phi, mu) {
    ### mu is a N x C matrix with ith row specifying mean of observation i
    ### row sums must be 1
    N <- nrow(mu)
    C <- ncol(mu)
    yout <- matrix(NA, nrow = N, ncol = C)
    for(i in 1:N) {
        yout[i,] <- gtools::rdirichlet(1, alpha = phi * mu[i,])
    }
    return(yout)
}

B1 <- matrix(c(.9, .05, .05,
               .05, .9, .05,
               .05, .09, .9), nrow = 3, byrow = TRUE)
B2 <- matrix(c(.4, .3, .3,
               .3, .4, .3,
               .3, .3, .4), nrow = 3, byrow = TRUE)
n <- 100
X <- rdirichlet(100, rep(1, 3))
mu <- X %*% B2
yout <- dir_data_gen(10, mu)


b_est <- pairedCompReg(yout = yout, ypred = X)
boot_out <- pairedCompReg(yout = yout, ypred = X, boot.ci = TRUE)
boot_out_mat <- boot_out$t
pp <- 0.95

ellipses_list <- lapply(1:3, function(i) {
    # find smallest ellipse containing the specified proportion of bootstrapped medians
    indices <- c(i, i + 3)
    fit <- cov.rob(boot_out_mat[,indices], quantile.used = ceiling(pp*nrow(boot_out_mat)), method = "mve")
    best_ellipse <- predict(ellipsoidhull(boot_out_mat[fit$best,indices] ))
    ellipse_df <- data.frame(y1 = best_ellipse[,1], y2 = best_ellipse[,2])
    ellipse_df$y3 <- pmax(0, 1 - ellipse_df$y1 -ellipse_df$y2)
    ellipse_df <- ellipse_df / rowSums(ellipse_df)
    return(ellipse_df)
})

est_df <- data.frame(y1 = b_est[,1],
                        y2 = b_est[,2],
                        y3 = b_est[,3])

### Plot 
ggtern(data = est_df , mapping = aes(x = y1, y = y2, z = y3)) +
    #geom_point()  +
    geom_mean_ellipse(data = ellipses_list[[1]]) +
    geom_mean_ellipse(data = ellipses_list[[2]]) +
    geom_mean_ellipse(data = ellipses_list[[3]])+
    geom_text(aes(label = lab)) +
    xlab(TeX("$E\\[y_1\\]$")) +
    ylab(TeX("$E\\[y_2\\]$")) +
    zlab(TeX("$E\\[y_3\\]$")) +
    ggtitle(TeX("$\\mathbf{B}^{1}")) +
    theme(plot.title = element_text(hjust = 0.5))
