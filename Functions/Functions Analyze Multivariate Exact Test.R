library(gridExtra)

#Statistical test
graphStatMultiExact=function(resSim_funct){
  
  sub_df=resSim_funct[resSim_funct$cm=="y ~ treatment + X1 + X2",]
  g1=ggplot(sub_df, aes(x = stat)) +
    geom_histogram(aes(y = ..density..), bins = 100, fill = "grey", color = "darkgrey") +
    stat_function(fun = stats::df, args = list(df1=1, df2=N-5, ncp=1), color = "black", size = 1) +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(
      title = expression(
        paste(
          "Test statistic"
        )
      ),
      x = NULL, 
      y = "Density" 
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),            
      axis.ticks.y = element_blank(),                    
      axis.text.y = element_blank()                    
    )+
    annotate(
      geom  = "text",
      x     = 8,
      y     = 0.12,
      label = "F[list(1, N-(K+1))](lambda==1)",  
      fontface = 2,
      parse = TRUE                    
    )
    #annotate(geom = "text",x=8,y=0.12,label=expression(paste('F'['1,N-(K+1)'],'(',lambda,'=1)')),fontface = 2)
  g1
  
  #Show that under the null the distribution of p-values is uniform:
  sub_df$pval=pf(sub_df$stat, df1=1, df2=N-5, ncp=1)
  g2=ggplot(sub_df, aes(x = pval)) +
    geom_histogram(aes(y = ..density..), breaks=seq(0,1,0.025), fill = "grey", color = "darkgrey")+
    stat_function(fun = dunif, args = list(min=0,max=1), color = "black", size = 1) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = expression(
        paste(
          #"Simulated distribution of p-values under the Null"
          "p-values"
        )
      ),
      x = NULL,
      y = NULL  
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size=10),            
      axis.title.y = element_blank(),                   
      axis.ticks.y = element_blank(),                    
      axis.text.y = element_blank()               
    )+
    coord_cartesian(ylim = c(0, 2))
  g=grid.arrange(g1, g2, ncol = 2)
  g
  
  ggsave(
    "Output/MultiExactCase_StatDistrib.pdf",
    g,
    width = 8,
    height = 4,
    dpi = 1200,
    bg='white'
  )
  
}

rejectionRateMultiExact=function(resSim_funct){
  
  sub_df=resSim_funct[resSim_funct$cm=="y ~ treatment + X1 + X2",]
  matRes_funct=matrix(data=NA,nrow=3,ncol=3)
  colnames(matRes_funct)=c("One-sided (low)", "One-sided (high)","Two-sided")
  rownames(matRes_funct)=c("1%","5%","10%")
  matRes_funct[1,1]=mean(sub_df$stat<=qf(0.01, df1=1, df2=N-5, ncp=1))
  matRes_funct[1,2]=mean(sub_df$stat>=qf(0.99, df1=1, df2=N-5, ncp=1))
  matRes_funct[1,3]=mean(sub_df$stat<=qf(0.005, df1=1, df2=N-5, ncp=1) | sub_df$stat>=qf(0.995, df1=1, df2=N-5, ncp=1))
  matRes_funct[2,1]=mean(sub_df$stat<=qf(0.05, df1=1, df2=N-5, ncp=1))
  matRes_funct[2,2]=mean(sub_df$stat>=qf(0.95, df1=1, df2=N-5, ncp=1))
  matRes_funct[2,3]=mean(sub_df$stat<=qf(0.025, df1=1, df2=N-5, ncp=1) | sub_df$stat>=qf(0.975, df1=1, df2=N-5, ncp=1))
  matRes_funct[3,1]=mean(sub_df$stat<=qf(0.10, df1=1, df2=N-5, ncp=1))
  matRes_funct[3,2]=mean(sub_df$stat>=qf(0.90, df1=1, df2=N-5, ncp=1))
  matRes_funct[3,3]=mean(sub_df$stat<=qf(0.05, df1=1, df2=N-5, ncp=1) | sub_df$stat>=qf(0.95, df1=1, df2=N-5, ncp=1))
  matRes_funct
}


