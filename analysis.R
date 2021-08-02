rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)
library(viridis)

setwd("~/Documents/GitHub/concensus_paper/")

data_files <- Sys.glob("Data_Files_Short/trial*csv")
column_names <- c("M_A", "M_R", "D", "M", "init_prop", "s", "radius", "delta", "mean_prop", "sd_prop", "min_prop", "max_prop", "prop_min", "prop_max", "kappa")

############################## DRAWING THE MOORE NEIGHBOURHOODS

moore_nbds <- rbind(
    reshape2::melt(matrix(rbind(
      c(0,0,0,0,0,0,0), 
      c(0,0,0,0,0,0,0), 
      c(0,0,0,0,0,0,0), 
      c(0,0,0,1,0,0,0), 
      c(0,0,0,0,0,0,0), 
      c(0,0,0,0,0,0,0),
      c(0,0,0,0,0,0,0)
    ), nrow=7)) |> dplyr::mutate(radius=1),
    reshape2::melt(matrix(rbind(
      c(0,0,0,0,0,0,0), 
      c(0,0,0,0,0,0,0), 
      c(0,0,1,1,1,0,0), 
      c(0,0,1,1,1,0,0), 
      c(0,0,1,1,1,0,0),
      c(0,0,0,0,0,0,0), 
      c(0,0,0,0,0,0,0)
    ), nrow=7)) |> dplyr::mutate(radius=2),
    reshape2::melt(matrix(rbind(
      c(0,0,0,0,0,0,0), 
      c(0,1,1,1,1,1,0), 
      c(0,1,1,1,1,1,0), 
      c(0,1,1,1,1,1,0), 
      c(0,1,1,1,1,1,0), 
      c(0,1,1,1,1,1,0),
      c(0,0,0,0,0,0,0)
    ), nrow=7)) |> dplyr::mutate(radius=3)
  )

neighbourhoods <- ggplot(moore_nbds, aes(x=Var1, y=Var2, group=radius, fill=value)) +
  geom_tile() +
  facet_grid(.~radius, labeller=label_both) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    aspect.ratio = 1,
    strip.text = element_blank()
  ) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_gradientn(colours = c("red","green"), values = c(0,1))
  ggsave(neighbourhoods, file="moore_neighbourhoods.png", width=15, height=5.5)
  system(sprintf("convert %s -trim %s", "moore_neighbourhoods.png", "moore_neighbourhoods.png"))

Finals <- data.table(file=as.character(), prop=as.numeric())
for(file_name in Sys.glob("Data_Files_Time/*.csv"))
{
  Results <- fread(file_name)
  Finals <- rbind(
    Finals,
    list(file = file_name, prop = Results[time==max(time), mean(prop_acceptors)])
  )
}

########################################################## TIME SERIES

# no clusters, so no M
Time_Series <- rbindlist(lapply(
    c(
      "Data_Files_Time/trial_p_0.3_s_3.0_r_6_delta_1.0.csv",
      "Data_Files_Time/trial_p_0.7_s_9.0_r_9_delta_9.0.csv",
      "Data_Files_Time/trial_p_0.4_s_8.0_r_6_delta_1.0.csv",
      "Data_Files_Time/trial_p_0.3_s_3.0_r_9_delta_1.0.csv",
      "Data_Files_Time/trial_p_0.3_s_9.0_r_5_delta_3.0.csv",
      "Data_Files_Time/trial_p_0.3_s_3.0_r_3_delta_1.0.csv",
      "Data_Files_Time/trial_p_0.3_s_9.0_r_10_delta_3.0.csv",
      "Data_Files_Time/trial_p_0.5_s_5.0_r_1_delta_1.0.csv"
    ),
    fread
  ))  |>
  dplyr::rename(D=linear_dimension, init_prop=initial_acceptor_proportion, radius=interaction_radius) %>%
  .[, .(accept_mean=mean(prop_acceptors), accept_sd=sd(prop_acceptors)), by=.(M_A, M_R, kappa, D, init_prop, s, radius, delta, time)]

# plot 5a
colour_vector <- viridis_pal()(8)
Figure_5a <- ggplot(data.table(), aes(x=time, y=accept_mean, ymin=pmax(0, accept_mean-accept_sd), ymax=pmin(1, accept_mean+accept_sd))) +
    geom_ribbon(data = Time_Series[init_prop==0.3 & s==3 & radius==6 & delta==1], aes(fill=colour_vector[1]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.3 & s==3 & radius==6 & delta==1], colour=colour_vector[1]) +
    geom_ribbon(data = Time_Series[init_prop==0.7 & s==9 & radius==9 & delta==9], aes(fill=colour_vector[2]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.7 & s==9 & radius==9 & delta==9], colour=colour_vector[2]) +
    geom_ribbon(data = Time_Series[init_prop==0.4 & s==8 & radius==6 & delta==1], aes(fill=colour_vector[3]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.4 & s==8 & radius==6 & delta==1], colour=colour_vector[3]) +
    geom_ribbon(data = Time_Series[init_prop==0.3 & s==3 & radius==9 & delta==1], aes(fill=colour_vector[4]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.3 & s==3 & radius==9 & delta==1], colour=colour_vector[4]) +
    geom_ribbon(data = Time_Series[init_prop==0.3 & s==9 & radius==5 & delta==3], aes(fill=colour_vector[5]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.3 & s==9 & radius==5 & delta==3], colour=colour_vector[5]) +
    geom_ribbon(data = Time_Series[init_prop==0.3 & s==3 & radius==3 & delta==1], aes(fill=colour_vector[6]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.3 & s==3 & radius==3 & delta==1], colour=colour_vector[6]) +
    geom_ribbon(data = Time_Series[init_prop==0.3 & s==9 & radius==10 & delta==3], aes(fill=colour_vector[7]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.3 & s==9 & radius==10 & delta==3], colour=colour_vector[7]) +
    geom_ribbon(data = Time_Series[init_prop==0.5 & s==5 & radius==1 & delta==1], aes(fill=colour_vector[8]), alpha=0.1) +
      geom_line(data = Time_Series[init_prop==0.5 & s==5 & radius==1 & delta==1], colour=colour_vector[8]) +
    labs(x=bquote("Time ("*tau*")"), y=expression(atop("Proportion of acceptors", "("*tau*"=10000)"))) + #, tags="A"
    theme_bw() +
      theme(
        legend.title.align=0.5,
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
      ) +
    scale_x_continuous(expand = c(0, 0), breaks=seq(0, 15000, by=1000)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(
      values = colour_vector,
      labels = c(
        "(0.3, 3, 6, 1)", "(0.7, 9, 9, 9)", "(0.4, 8, 6, 1)", "(0.3, 3, 9, 1)", 
        "(0.3, 9, 5, 3)", "(0.3, 3, 3, 1)", "(0.3, 9, 10, 3)", "(0.5, 5, 1, 1)"
      )
    ) +
    guides(fill=guide_legend(title=bquote("(p, s, r, "*delta*")")))
    ggsave(Figure_5a, file="spatial_time_series.png", width=10, height=3, dpi=600)
    
    
########################################### RESULTS GRAPHS
    
Finals_Mean_and_SD <- data.table(matrix(nrow=length(data_files), ncol=length(column_names), -1.0))
names(Finals_Mean_and_SD) <- column_names

start_time <-  Sys.time()
for(i in seq(data_files))
{
    Results <- fread(data_files[i]) |>
      dplyr::rename(D=network_linear_dimension, M=cluster_linear_dimension, init_prop=initial_acceptor_proportion, radius=interaction_radius)

    Results <- Results[
            Results[, .I[time==max(time)], by=trial]$V1
        ][,
          .(
            mean_prop=mean(prop_acceptors), sd_prop=sd(prop_acceptors),
            min_prop=min(prop_acceptors), max_prop=max(prop_acceptors)
          ),
          by=.(M_A, M_R, D, M, init_prop, s, radius, delta)
        ] |>
        dplyr::mutate(
          prop_min = pmax(0, mean_prop-sd_prop),
          prop_max = pmin(1, mean_prop+sd_prop),
          kappa:=M_R/M_A
      )

    Finals_Mean_and_SD[i, (colnames(Finals_Mean_and_SD)):=as.list(Results)]
}
Finals_Mean_and_SD <- Finals_Mean_and_SD |>
  dplyr::mutate(
    M_A = factor(M_A),
    M_R = factor(M_R),
    dim = factor(dim),
    init_prop = factor(init_prop),
    delta = factor(delta),
  )
print(Sys.time() - start_time)
fwrite(Finals_Mean_and_SD,  "total_t_10000_means_and_sds.csv")

Finals_Mean_and_SD <- fread("total_t_10000_means_and_sds.csv")
    
Figure_5b <- ggplot(
        Finals_Mean_and_SD[delta %in% c(0, 0.5, 1) & s==10 & M==0] |> 
            # dplyr::mutate(delta = unlist(lapply(delta, \(x) expression(delta*"= ") )) ) +
            dplyr::mutate(init_prop=factor(init_prop), delta=factor(delta)),
        aes(x=radius, y=mean_prop, group=init_prop, colour=init_prop)
    ) +
    geom_line(size=1) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=pmax(0, mean_prop-sd_prop), ymax=pmin(1, mean_prop+sd_prop))) +
    theme_bw() +
    theme(
        legend.title.align = 0.5,
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
    ) +
    labs(x=bquote("Interaction radius ("*r*")"), y=expression(atop("Proportion of acceptors", "("*tau*"=10000)")) ) +
    guides(colour=guide_legend(title="p")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_colour_viridis_d() +
    facet_grid(.~delta, labeller=as_labeller(c( "0"="delta*' = 0.0'", "0.5"="delta*' = 0.5'", "1"="delta*' = 1.0'" ), default=label_parsed)) +
    ylim(c(0,1))
    ggsave(Figure_5b, file="spatial_deltas.png", height=3, width=10)

Figure_5d <- ggplot(
        Finals_Mean_and_SD[kappa==0.34 & s%in%0:10 & delta%in%seq(0,5,by=0.5) & radius%in%1:25 & M==0] |> 
          dplyr::mutate(s=sprintf("s %02i", s), delta=sprintf("d %.1f", delta)),
        aes(x=radius, y=mean_prop, group=init_prop, colour=init_prop)
    ) +
    geom_line(size=1) +
    scale_colour_viridis() +
    geom_ribbon(aes(ymin=min_prop, ymax=max_prop), alpha=0.1, linetype=1, colour="grey70") +
    facet_grid(
      s~delta,
      labeller=as_labeller(c(
        "d 0.0" = "delta*' = 0.0'",
        "d 0.5" = "delta*' = 0.5'",
        "d 1.0" = "delta*' = 1.0'",
        "d 1.5" = "delta*' = 1.5'",
        "d 2.0" = "delta*' = 2.0'",
        "d 2.5" = "delta*' = 2.5'",
        "d 3.0" = "delta*' = 3.0'",
        "d 3.5" = "delta*' = 3.5'",
        "d 4.0" = "delta*' = 4.0'",
        "d 4.5" = "delta*' = 4.5'",
        "d 5.0" = "delta*' = 5.0'",
        "s 00" = "s*' = 0'",
        "s 01" = "s*' = 1'",
        "s 02" = "s*' = 2'",
        "s 03" = "s*' = 3'",
        "s 04" = "s*' = 4'",
        "s 05" = "s*' = 5'",
        "s 06" = "s*' = 6'",
        "s 07" = "s*' = 7'",
        "s 08" = "s*' = 8'",
        "s 09" = "s*' = 9'",
        "s 10" = "s*' = 10'"
      ), default=label_parsed)
    ) +
    theme_bw() +
    labs(
      x="Interaction radius (r)", 
      y=bquote("Proportion of acceptors ("*tau*"=10000)")
    ) +
    theme(
        legend.title.align=0.5,
        strip.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        panel.spacing.y=unit(1, "lines")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    guides(colour=guide_legend(bquote(p)))
  ggsave(Figure_5d, file="spatial_appendix.png", height=20, dpi=600)
