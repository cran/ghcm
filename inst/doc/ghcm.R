## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  dpi=125,
  fig.width=6,
  fig.height=4,
  out.width="100%"
)
options(rmarkdown.html_vignette.check_title = FALSE) 

## ----setup, include=FALSE-----------------------------------------------------
library(ghcm)

## -----------------------------------------------------------------------------
library(ghcm)
set.seed(111)
data(ghcm_sim_data)
grid <- seq(0, 1, length.out=101)
colnames(ghcm_sim_data)

## ----plot-X, fig.cap="Plot of $X$ with the estimated mean curve in red.", echo=FALSE----
matplot(grid, t(ghcm_sim_data$X), type="l",
        col=rgb(0,0,0,0.1) , lty=1, xlab="",
        ylab="")
lines(grid, colMeans(ghcm_sim_data$X), col=2,
      lty=2, lwd=2)
title("X")

## ----plot-Z, fig.cap="Plot of $Z$ with the estimated mean curve in red.", echo=FALSE----
matplot(grid, t(ghcm_sim_data$Z), type="l",
        col=rgb(0,0,0,0.1) , lty=1, xlab="",
        ylab="")
lines(grid, colMeans(ghcm_sim_data$Z), col=2,
      lty=2, lwd=2)
title("Z")

## ----plot-W, fig.cap="Plot of $W$ with the estimated mean curve in red.", echo=FALSE----
matplot(grid, t(ghcm_sim_data$W), type="l",
        col=rgb(0,0,0,0.1) , lty=1, xlab="",
        ylab="")
lines(grid, colMeans(ghcm_sim_data$W), col=2,
      lty=2, lwd=2)
title("W")

## -----------------------------------------------------------------------------
cor(ghcm_sim_data$Y_1, ghcm_sim_data$Y_2)

## ----scalar-given-functional-test---------------------------------------------
library(refund)
m_1 <- pfr(Y_1 ~ lf(X) + lf(Z) + lf(W) , data=ghcm_sim_data)
m_2 <- pfr(Y_2 ~ lf(X) + lf(Z) + lf(W), data=ghcm_sim_data)
test <- ghcm_test(resid(m_1), resid(m_2))
print(test)

## ----Z-Y-plot, echo=FALSE, fig.cap="Plot of $Z$ with colors based on the value of $Y_1$."----
library(ggplot2)
library(reshape2)

Z_df <- as.data.frame(ghcm_sim_data$Z)
colnames(Z_df) <- grid
Z_df$Y_1<- ghcm_sim_data$Y_1
Z_df$id <- 1:nrow(ghcm_sim_data)
tmp <- melt(Z_df, id.var=c("id", "Y_1"), variable.name="grid")
tmp$grid <- as.numeric(levels(tmp$grid))[tmp$grid]


ggplot(tmp, aes(x=grid, y=value, group=id)) + geom_line(aes(color=Y_1)) +
  theme_bw() + scale_x_continuous(name="") + scale_y_continuous(name="") +
  scale_color_gradient(name=expression(Y[1]))

## ----X-Y-plot, echo=FALSE, fig.cap="Plot of $X$ with colors based on the value of $Y_1$."----
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

X_df <- as.data.frame(ghcm_sim_data$X)
colnames(X_df) <- grid
X_df$Y_1<- ghcm_sim_data$Y_1
X_df$id <- 1:nrow(ghcm_sim_data)
tmp <- melt(X_df, id.var=c("id", "Y_1"), variable.name="grid")
tmp$grid <- as.numeric(levels(tmp$grid))[tmp$grid]


ggplot(tmp, aes(x=grid, y=value, group=id)) + geom_line(aes(color=Y_1)) +
  theme_bw() + scale_x_continuous(name="") + scale_y_continuous(name="") +
  scale_color_gradient(name=expression(Y[1]))

## ----Y-X-Z-test---------------------------------------------------------------
m_1 <- pfr(Y_1 ~ lf(Z), data = ghcm_sim_data)
m_X <- pffr(X ~ ff(Z), data = ghcm_sim_data, chunk.size = 31000)
test <- ghcm_test(resid(m_X), resid(m_1))
print(test)

## ----X-W-Z-test---------------------------------------------------------------
m_X <- pffr(X ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
m_W <- pffr(W ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
test <- ghcm_test(resid(m_X), resid(m_W))
print(test)

## -----------------------------------------------------------------------------
data(ghcm_sim_data_irregular)

## -----------------------------------------------------------------------------
head(ghcm_sim_data_irregular$X)

## ----plot-X-irregular, fig.cap="Plot of the first 5 observations of $X$ with the subsampled grid points marked by the larger triangles.", echo=FALSE----
regular_X <- melt(ghcm_sim_data$X, varnames = c(".obs", ".index"), value.name = ".value")
regular_X$.index <- grid[regular_X$.index]
plot_df <- left_join(regular_X,ghcm_sim_data_irregular$X %>% mutate(subsampled=TRUE)) %>% replace_na(list(subsampled=FALSE))
plot_df %>% subset(.obs <= 5) %>% ggplot(aes(x=.index, y=.value, color=as.factor(.obs))) + geom_point(aes(shape=subsampled, size=subsampled)) + geom_line(alpha=0.3) +
  scale_size_manual(values=c(0.5, 2)) + guides(color="none", size="none", shape="none") +
  scale_x_continuous(name=NULL) + scale_y_continuous(name=NULL) + theme_bw()

## ----plot-W-irregular, fig.cap="Plot of the first 5 observations of $W$ with the subsampled grid points marked by the larger triangles.", echo=FALSE----
regular_W <- melt(ghcm_sim_data$W, varnames = c(".obs", ".index"), value.name = ".value")
regular_W$.index <- grid[regular_W$.index]
plot_df <- left_join(regular_W,ghcm_sim_data_irregular$W %>% mutate(subsampled=TRUE)) %>% replace_na(list(subsampled=FALSE))
plot_df %>% subset(.obs <= 5) %>% ggplot(aes(x=.index, y=.value, color=as.factor(.obs))) + geom_point(aes(shape=subsampled, size=subsampled)) + geom_line(alpha=0.3) +
  scale_size_manual(values=c(0.5, 2)) + guides(color="none", size="none", shape="none") +
  scale_x_continuous(name=NULL) + scale_y_continuous(name=NULL) + theme_bw()

## ----Y-X-Z-test-irregular-----------------------------------------------------
n <- nrow(ghcm_sim_data_irregular$Z)
Z_df <- data.frame(.obs=1:n)
Z_df$Z <- ghcm_sim_data_irregular$Z
m_1 <- pfr(Y_1 ~ lf(Z), data = ghcm_sim_data_irregular)
m_X <- pffr(X ~ ff(Z), ydata = ghcm_sim_data_irregular$X, data=Z_df, chunk.size=31000)
test <- ghcm_test(resid(m_X), resid(m_1), X_limits=c(0, 1))
print(test)

## ----X-W-Z-test-irregular-----------------------------------------------------
n <- nrow(ghcm_sim_data_irregular$Z)
Z_df <- data.frame(.obs=1:n)
Z_df$Z <- ghcm_sim_data_irregular$Z

m_X <- pffr(X ~ ff(Z), ydata = ghcm_sim_data_irregular$X, data=Z_df, chunk.size=31000)
m_W <- pffr(W ~ ff(Z), ydata = ghcm_sim_data_irregular$W, data=Z_df, chunk.size=31000)
test <- ghcm_test(resid(m_X), resid(m_W), X_limits=c(0, 1), Y_limits=c(0, 1))
print(test)

