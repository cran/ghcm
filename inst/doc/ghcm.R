## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  dpi=150,
  fig.width=7,
  fig.height=5,
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
test <- ghcm_test(resid(m_1), resid(m_2), X_grid = NA, Y_grid = NA )
print(test)

## ----scalar-given-functional-test-plot, echo=FALSE, fig.cap="Plot of the estimated asymptotic test distribution under the null with the red line indicating the observed value."----
plot(test)

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
test <- ghcm_test(resid(m_X), resid(m_1), X_grid = grid, Y_grid = NA )
print(test)

## ----Y-X-Z-test-plot, echo=FALSE, fig.cap="Plot of the estimated asymptotic test distribution under the null with the red line indicating the observed value."----
plot(test)

## ----X-W-Z-test---------------------------------------------------------------
m_X <- pffr(X ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
m_W <- pffr(W ~ ff(Z), data=ghcm_sim_data, chunk.size=31000)
test <- ghcm_test(resid(m_X), resid(m_W), X_grid = grid, Y_grid = grid )
print(test)

## ----X-W-Z-test-plot, echo=FALSE, fig.cap="Plot of the estimated asymptotic test distribution under the null with the red line indicating the observed value."----
plot(test)

