library(ggplot2)

#Statistical test
graphStatUni=function(resSim_funct){
  
  resSim_funct$stat=(resSim_funct$gamma/resSim_funct$gamma_se)^2
  g1=ggplot(resSim_funct, aes(x = stat)) +
    geom_histogram(aes(y = ..density..), bins = 50, fill = "grey", color = "darkgrey") +
    stat_function(fun = stats::df, args = list(df1=1, df2=N-2, ncp=1), color = "black", size = 1) +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(
      title = expression(
        paste(
          # "Simulated distribution of test statistics (", 
          # t[hat(gamma)], 
          # ") under the Null"
          "Test statistic"
        )
      ),
      x = NULL, # Leave x-label empty
      y = "Density"  # Leave y-label empty
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),            # Center the title
      axis.ticks.y = element_blank(),                    # Remove y-axis ticks
      axis.text.y = element_blank()                     # Remove y-axis text
    )+
    annotate(
      geom  = "text",
      x     = 6,
      y     = 0.12,
      label = "F[list(1, N-2)](lambda==1)",
      fontface = 2,
      parse = TRUE                    
    )
  g1
  
  resSim_funct$pval=pf(resSim_funct$stat, df1=1,df2=N-2,ncp=1)
  g2=ggplot(resSim_funct, aes(x = pval)) +
    geom_histogram(aes(y = ..density..), breaks=seq(0,1,0.05), fill = "grey", color = "darkgrey")+
    stat_function(fun = dunif, args = list(min=0,max=1), color = "black", size = 1) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = expression(
        paste(
          "p-values"
        )
      ),
      x = NULL, # Leave x-label empty
      y = NULL  # Leave y-label empty
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size=10),            # Center the title
      axis.title.y = element_blank(),                    # Remove y-axis label
      axis.ticks.y = element_blank(),                    # Remove y-axis ticks
      axis.text.y = element_blank()                     # Remove y-axis text
    )+
    coord_cartesian(ylim = c(0, 2))
  g=grid.arrange(g1, g2, ncol = 2)
  g
  
  ggsave(
    "Output/UnivCase_StatDistrib.pdf",
    g,
    width = 8,
    height = 4,
    dpi = 1200,
    bg='white'
  )
  
}

rejectionRateUni=function(resSim_funct){
  resSim_funct$stat=(resSim_funct$gamma/resSim_funct$gamma_se)^2
  matRes_funct=matrix(data=NA,nrow=3,ncol=3)
  colnames(matRes_funct)=c("One-sided (low)", "One-sided (high)","Two-sided")
  rownames(matRes_funct)=c("1%","5%","10%")
  matRes_funct[1,1]=mean(resSim_funct$stat<=qf(0.01, df1=1, df2=N-2, ncp=1))
  matRes_funct[1,2]=mean(resSim_funct$stat>=qf(0.99, df1=1, df2=N-2, ncp=1))
  matRes_funct[1,3]=mean(resSim_funct$stat<=qf(0.005, df1=1, df2=N-2, ncp=1) | resSim_funct$stat>=qf(0.995, df1=1, df2=N-2, ncp=1))
  matRes_funct[2,1]=mean(resSim_funct$stat<=qf(0.05, df1=1, df2=N-2, ncp=1))
  matRes_funct[2,2]=mean(resSim_funct$stat>=qf(0.95, df1=1, df2=N-2, ncp=1))
  matRes_funct[2,3]=mean(resSim_funct$stat<=qf(0.025, df1=1, df2=N-2, ncp=1) | resSim_funct$stat>=qf(0.975, df1=1, df2=N-2, ncp=1))
  matRes_funct[3,1]=mean(resSim_funct$stat<=qf(0.10, df1=1, df2=N-2, ncp=1))
  matRes_funct[3,2]=mean(resSim_funct$stat>=qf(0.90, df1=1, df2=N-2, ncp=1))
  matRes_funct[3,3]=mean(resSim_funct$stat<=qf(0.05, df1=1, df2=N-2, ncp=1) | resSim_funct$stat>=qf(0.95, df1=1, df2=N-2, ncp=1))
  matRes_funct
}


