# Code to produce surface plots of Figure 1.

# Define relevant functions
expit <- function(x) exp(x)/(1+exp(x))
m <- function(z1,z2,z3) (expit(z1)+expit(z2)+expit(z3)) * (0.9-0.8*expit(z1))
h <- function(z1,z2,z3) {
    mu = expit(z1)+expit(z2)+expit(z3)
    V = 0.1
    k = mu^2/V
    theta = V/mu
    alpha = 1/0.4*pgamma(1.5, shape=k+1, scale=theta) + 1/0.1*(1-pgamma(1.5, shape=k+1, scale=theta))
    alpha_test = mu * alpha
    beta_test = 1/0.4*pgamma(1.5, shape=k, scale=theta) + 1/0.1*(1-pgamma(1.5, shape=k, scale=theta))
    P = 0.1+0.8*expit(z1)
    hh = (alpha*(1-P)) / (beta_test*(1-P)+1/0.9*P)
    return(hh)
}

Z1 <- -1
grid <- expand.grid(x = seq(-10,10,0.1), y = seq(-10,10,0.1))
M <- function(x,y) m(Z1,x,y)
G <- function(x,y) h(Z1,x,y)


library(ggplot2)
library(gridExtra)
library(latex2exp)
plot_m = ggplot(grid, aes(x, y, z = M(x,y))) +
  geom_contour_filled(color = "black") +
  labs(x=bquote(Z[2]), y=bquote(Z[3]), title=expression(paste(m[0](Z)==bgroup("", bold(E) * '[' * 'X' * ' | ' * Z * ']', "")))) +
  labs(colour = "") +
  theme(axis.title=element_text(size=14))

plot_h = ggplot(grid, aes(x, y, z = G(x,y))) +
  geom_contour_filled(color = "black") +
  labs(x=bquote(Z[2]), y=bquote(Z[3]), title=expression(paste(h[0](Z)== bgroup("", bold(E) * '[' * sigma[0]^{-2} * '(X,Z)' * ' | ' * Z * ']', "")^{-1}, bgroup("", bold(E) * '[' * sigma[0]^{-2} * '(X,Z) X' * ' | ' * Z * ']', "")))) +
  labs(colour = "") +
  theme(axis.title=element_text(size=14))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
grid.arrange(plot_m+theme(legend.position="none"), plot_h+theme(legend.position="none"), get_legend(plot_m), ncol=3, widths=c(3, 3, 1.2))

# For contour plots (note each contour plot here has a different colour scaling!)
grid.arrange(plot_m+theme(legend.position="none"), plot_h+theme(legend.position="none"), get_legend(plot_m), ncol=3, widths=c(3, 3, 1.2))



# Code for surface plots:

x <- seq(-10, 10, length = 30)
y <- seq(-10, 10, length = 30)
z <- outer(x, y, function(x, y) G(x,y))
z_matrix <- matrix(z, nrow = length(x), ncol = length(y))
plot_ly(x = ~x, y = ~y, z = ~z_matrix) %>%
  add_surface() %>%
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = 1))))


library(plotly)
library(htmltools)
x1 <- seq(-4, 4, length = 300)
y1 <- seq(-4, 4, length = 300)
z1 <- outer(x1, y1, function(x, y) M(x,y))
z_matrix1 <- matrix(z1, nrow = length(x1), ncol = length(y1))
x2 <- seq(-4, 4, length = 300)
y2 <- seq(-4, 4, length = 300)
z2 <- outer(x2, y2, function(x, y) G(x,y))
z_matrix2 <- matrix(z2, nrow = length(x2), ncol = length(y2))
# m(Z) surface
plot1 <- plot_ly(x = ~x1, y = ~y1, z = ~z_matrix1, showscale = FALSE) %>%
  add_surface() %>%
  layout(scene = list(
    xaxis = list(title = TeX("$Z_1$")),
    yaxis = list(title = "Y Axis"),
    zaxis = list(title = "Z Axis"),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1))
  ),
  title = "m")
# h(Z) surface
plot2 <- plot_ly(x = ~x2, y = ~y2, z = ~z_matrix2, showscale = FALSE) %>%
  add_surface() %>%
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = 1))),
         title = bquote(h[0](X)))
# Both plots (note plot is interactive)
browsable(
  tagList(
    div(style = "display: inline-block; width: 45%;", plot1),
    div(style = "display: inline-block; width: 45%;", plot2)
  )
)
