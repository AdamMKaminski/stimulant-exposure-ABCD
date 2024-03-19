
# themods: a list of bayesian models to extract posteriors from
# range: desired predictors in model
# titl: title for plot

plot_RandEffPosteriors <- function(themods,range,titl) {
  
  iter <- 1
  for (mod in themods) {
    temp <- posterior_summary(mod,prob=c(0.005, .025, .975, 0.995))
    temp2 <- as.data.frame(temp)
    temp3 <-as.data.frame(temp2[range,])
    temp4 <- as.data.frame(t(temp3))
    temp3$fc <- colnames(temp4)
    if (iter==1){
      tograph <- temp3
    } else {
      tograph <- rbind(tograph,temp3)
    }
    iter <- iter + 1
  }
  
  #
  # FOR 10 NETWORKS:
  #
  
  tograph$network <- c("L caudate","L caudate","L caudate","L caudate","L caudate","L caudate","L caudate","L caudate","L caudate","L caudate",
                       "L putamen","L putamen","L putamen","L putamen","L putamen","L putamen","L putamen","L putamen","L putamen","L putamen",
                       "L NAcc","L NAcc","L NAcc","L NAcc","L NAcc","L NAcc","L NAcc","L NAcc","L NAcc","L NAcc",
                       "R caudate","R caudate","R caudate","R caudate","R caudate","R caudate","R caudate","R caudate","R caudate","R caudate",
                       "R putamen","R putamen","R putamen","R putamen","R putamen","R putamen","R putamen","R putamen","R putamen","R putamen",
                       "R NAcc","R NAcc","R NAcc","R NAcc","R NAcc","R NAcc","R NAcc","R NAcc","R NAcc","R NAcc")
  
  require(forcats)
  require(ggplot2)
  require(cowplot)
  
  tograph$sig <- as.numeric((tograph$Q2.5>0&tograph$Q97.5>0)|(tograph$Q2.5<0&tograph$Q97.5<0))
  a <- ifelse((tograph$sig==1), "red","black")
  
  ggplot(tograph,aes(x = as.factor(fc), y = Estimate, color = as.factor(network))) +
    theme_cowplot(12) +
    geom_point(position=position_dodge(width=0.5)) +
    aes(x = fct_inorder(as.factor(fc))) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = I(0.2), color = I("black")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color=a)) +
    geom_errorbar(aes(x = as.factor(fc),
                      ymin = Q0.5,
                      ymax = Q99.5),
                  width = .1,
                  position=position_dodge(width=0.5)) +
    labs(color = "Seeds", x="Networks", y="Posterior Distributions for Striatal rs-FC") +
    theme(legend.position = c(0.45, 0.1),
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold")
    ) +
    scale_color_discrete(breaks=c('L caudate', 'L putamen', 'L NAcc', 'R caudate', 'R putamen', 'R NAcc')) +
    geom_point(shape='-', size=5, color = 'black', data = data.frame(x=tograph$fc,
                                                                     y=tograph$Q2.5,
                                                                     color=tograph$network), aes(x=x, y=y, color=color)) +
    geom_point(shape='-', size=5, color = 'black', data = data.frame(x=tograph$fc,
                                                                     y=tograph$Q97.5,
                                                                     color=tograph$network), aes(x=x, y=y, color=color)) +
    scale_x_discrete(labels=rep(c("Auditory","Cingulo-op-salience","Cingulo-parietal","Default Mode","Dorsal Attention","Fronto-parietal","Retrosplenial","Somatomotor","Ventral Attention","Visual"),times=6)) + 
    ggtitle(titl) +
    theme(plot.title = element_text(size = 12.5))
  
}