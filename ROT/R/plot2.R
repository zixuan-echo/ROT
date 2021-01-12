require(gridExtra)
require(ggplot2)
#'@title Plot 1D distribution
#'@description Plot matrix M with the source and target 1D distribution
#'Creates a subplot with the source distribution a on the left and target distribution b on the tot. The matrix M is shown in between.
#'
#'@param a Source distribution
#'@param b Target distribution
#'@param M Matrix to plot
#'@param title The title of the plot
#'
#' @examples 
#' ## SMALL EXAMPLE
#' n=10
#' a=dnorm(1:100,20,10)
#' a=a/sum(a)
#' b=dnorm(1:100, mean=60, sd=2)
#' b=b/sum(b)*2
#' distance= dists(1:100,1:100)
#' distance= diatance/max(distance)
#' plan=sinkhorn_knopp_unbalanced(a,b,distance,1,1,1000,1e-6,TRUE,FALSE)
#' plot1D(a, b, plan, title = "UOT example")
plot1D <- function(a, b, M, title = ""){
  na = length(a)
  nb = length(b)
  ap = apply(M,1,sum)
  bp = apply(M,2,sum)
  dat_a = data.frame(ind = 1:na, a, ap)
  dat_b = data.frame(ind = 1:nb, b, bp)
  hist_top <- ggplot(dat_b) +
    geom_line(aes(x = ind, y = b), colour="red",size=1)+
    geom_line(aes(x = ind, y = bp), colour="pink",size=1)+
    ggtitle(title)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0))+
    theme(plot.title=element_text(size=15,face="plain",hjust=.5),
          panel.background=element_blank(),
          axis.title = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank())
  hist_right <- ggplot(dat_a) +
    geom_line(aes(x = ind, y = a), colour="blue",size=1)+
    geom_line(aes(x = ind, y = ap), colour="lightblue",size=1)+
    scale_x_reverse(expand = c(0,0))+
    scale_y_continuous(position = "left", expand = c(0,0))+
    coord_flip()+
    theme(panel.background=element_blank(),
          axis.title = element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank())
  row_ind = rep(dim(M)[1]:1, each = dim(M)[2], times = 1)
  col_ind = rep(1:dim(M)[2], each = 1, times = dim(M)[1])
  m = as.vector(t(M))
  dat_m = data.frame(row_ind, col_ind, m)
  scatter <- ggplot(dat_m,aes(x=col_ind,y=row_ind,fill=m))+
    geom_raster()+ scale_fill_gradient2(low="red", high="blue", mid="white")+
    guides(fill=F)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0))+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
          axis.title  = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.line = element_blank()
    )
  empty <- ggplot()+theme_minimal()
  grid.arrange(empty, hist_top, hist_right, scatter, ncol=2, nrow=2, widths=c(1.2,4), heights=c(1.5,4))
}

#'@title Plot 2D distribution
#'@description Plot matrix M in 2D with lines using alpha values.
#'@param a an \eqn{(M\times 2)} matrix of Source samples positions
#'@param b an \eqn{(N\times 2)} matrix of Target samples positions
#'@param M OT matrix
#'@param thr threshold above which the line is drawn
#'@param title The title of the plot
#'
#' @examples 
#' ## SMALL EXAMPLE
#' n = 30
#' mu_s = c(0, 0)
#' cov_s = matrix(c(1, 0, 0, 1), 2, 2)
#' mu_t = c(4, 4)
#' cov_t = matrix(c(1, -0.5, -0.5, 1), 2, 2)
#' library(MASS)
#' library(transport)
#' xs = mvrnorm(n, mu_s, cov_s)
#' xt = mvrnorm(n, mu_t, cov_t)
#' a = rep(1, n)/n
#' b = rep(1, n)/n
#' M = dist_mat(xs, xt)
#' M = M/max(M)
#' plan = transport(a, b, M)
#' p = rep(0, n^2)
#' p[(plan$to-1)*n + plan$from] = plan$mass
#' p = matrix(p, n, n)
#' plot2D(xs, xt, p, title = "2DOT example")
plot2D <- function(xs, xt, G, thr=1e-8, title = ""){
  t = data.frame(
    x = numeric(0),
    y = numeric(0),
    group = numeric(0),
    alpha = numeric(0)
  )
  ds = as.data.frame(xs)
  dt = as.data.frame(xt)
  for(i in 1:dim(G)[1]){
    for(j in 1:dim(G)[2]){
      if(G[i, j]>thr){
        t <- rbind(t, c(xs[i, 1], xs[i, 2], g, G[i, j]), c(xt[j, 1], xt[j, 2], g, G[i, j]))
        g <- g + 1
      }
    }
  }
  colnames(t) <- c("x", "y", "group", "alpha")
  t$group <- as.factor(t$group)
  ggplot()+geom_line(data = t, aes(x = x, y = y, group = group, alpha = alpha), color = rgb(0.5 ,0.5 , 1))+
    geom_point(data = ds, aes(x = V1, y = V2), shape = 3, color = "blue")+
    geom_point(data = dt, aes(x = V1, y = V2), shape = 4, color = "red")+
    ggtitle(title)+
    guides(alpha=F)+
    theme(
      plot.title=element_text(size=15,face="plain",hjust=.5),
      panel.background=element_rect(fill="white",colour="black",size=0.25),
      axis.line=element_line(colour="black",size=0.25),
      axis.title=element_text(size=13,face="plain",color="black"),
      axis.title.x  = element_blank(),
      axis.title.y  = element_blank(),
    )
}