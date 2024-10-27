# Packages ----
# Ensures the package "pacman" is installed
if (!require("pacman")) install.packages("pacman")

pacman::p_load(dplyr,
               metafor,
               bayesmeta,
               brms,
               here,
               tidybayes,
               tidyr,
               ggplot2,
               MetBrewer,
               patchwork)

renv::restore()

# Data ----

d = dplyr::tribble(
  ~outcome, ~type, ~study, ~events_short, ~total_short, ~events_long, ~total_long,
  "mortality30", "RCT", "Molina, 2022", 3, 119, 9, 129,
  "mortality30", "RCT", "von Dach, 2020", 6, 169, 4, 165,
  "mortality30", "RCT", "Yahav, 2019", 15, 306, 13, 298,
  
  "mortality30", "obs", "Bae, 2021", 6, 97, 18, 193,
  "mortality30", "obs", "Chotiprasitsakul, 2018", 37, 385, 39, 385,
  "mortality30", "obs", "Fabre, 2019", 5, 72, 6, 179,
  "mortality30", "obs", "Giannella, 2018", 11, 426, 13, 430,
  "mortality30", "obs", "Ruiz-Ruigomez, 2020", 3, 23, 1, 29,
  "mortality30", "obs", "Sousa, 2019", 23, 163, 23, 232,
  "mortality30", "obs", "Surapat, 2020", 16, 41, 19, 110,
  
  "relapse30", "RCT", "Molina, 2022", 7, 119, 6, 129,
  "relapse30", "RCT", "Von Dach, 2020", 1, 169, 2, 165,
  
  "relapse30", "obs", "Bae, 2021", 5, 97, 15, 193,
  "relapse30", "obs", "Chotiprasitsakul, 2018", 5, 385, 9, 385,
  "relapse30", "obs", "Fabre, 2019", 5, 72, 20, 179,
  "relapse30", "obs", "Ruiz-Ruigomez, 2020", 1, 23, 0, 29,
  "relapse30", "obs", "Surapat, 2020", 2, 41, 7, 110
)

d_logOR = 
  metafor::escalc(
    measure = "OR", # log odds ratio,
    
    # short course
    ai = events_short,
    n1i = total_short,
    
    # Long course
    ci = events_long,
    n2i = total_long,
    
    data = d
  ) |> 
  dplyr::mutate(sei = sqrt(vi),
                type = factor(type, levels = c("RCT", "obs")))

# Models ----

mf = 
  formula(yi | se(sei) ~ 0 + type + (1|study))

# Tau
informative = bayesmeta::TurnerEtAlPrior("all-cause mortality",
                                         "pharmacological",
                                         "placebo / control")


logmean = informative$parameters["tau", "mu"]
logsd = informative$parameters["tau", "sigma"]

# logmean
# -1.975

# logsd
# 0.67


priors_mortality = 
  brms::prior(normal(0, 0.82), class = "b") +
  brms::prior(lognormal(-1.975, 0.67), class = "sd") 

m_mortality = 
  brm(
    data = subset(d_logOR, outcome == "mortality30"),
    family = gaussian,
    
    formula = mf,
    prior = priors_mortality,
    
    control = list(adapt_delta = .97),
    backend = "cmdstanr", # faster
    cores = parallel::detectCores(),
    chains = 4,
    warmup = 2000,
    iter = 4000,
    seed = 123,
    
    file = here::here("fits", "m_mortality.Rds"),
    file_refit = "on_change"
  )

# Tau
informative = bayesmeta::TurnerEtAlPrior(
  "signs / symptoms reflecting continuation / end of condition",
  "pharmacological",
  "placebo / control")

logmean = informative$parameters["tau", "mu"]
logsd = informative$parameters["tau", "sigma"]

# logmean
# -1.03

# logsd
# 0.755

priors_relapse = 
  brms::prior(normal(0, 0.82), class = "b") +
  brms::prior(lognormal(-1.03, 0.755), class = "sd") 

m_relapse = 
  brm(
    data = subset(d_logOR, outcome == "relapse30"),
    family = gaussian,
    
    formula = mf,
    prior = priors_relapse,
    
    control = list(adapt_delta = .97),
    backend = "cmdstanr", # faster
    cores = parallel::detectCores(),
    chains = 4,
    warmup = 2000,
    iter = 4000,
    seed = 123,
    
    
    file = here::here("fits", "m_relapse.Rds"),
    file_refit = "on_change"
  )

# Results wrangling ----

mortality_draws = 
  m_mortality |> 
  tidybayes::tidy_draws() |>
  dplyr::summarise(tau = sd_study__Intercept,
                   RCT = b_typeRCT,
                   obs = b_typeobs,
                   logROR = RCT - obs) |> 
  tidyr::pivot_longer(1:4) 

mortality_hdi = 
  mortality_draws |> 
  dplyr::group_by(name) |> 
  tidybayes::median_hdi()

relapse_draws = 
  m_relapse |> 
  tidybayes::tidy_draws() |>
  dplyr::summarise(tau = sd_study__Intercept,
                   RCT = b_typeRCT,
                   obs = b_typeobs,
                   logROR = RCT - obs) |> 
  tidyr::pivot_longer(1:4) 

relapse_hdi = 
  relapse_draws |> 
  dplyr::group_by(name) |> 
  tidybayes::median_hdi()  

# Forest plots ----  

pal = MetBrewer::met.brewer(name="VanGogh1",n=7,type="discrete")

forest_mortality = function(){
  
  dd = subset(d_logOR, outcome == "mortality30") |> 
    dplyr::arrange(desc(yi))
  
  x = 
    base::with(
      dd,
      
      metafor::forest(x = yi,
                      sei = sei,
                      slab = study,
                      order = desc(type),
                      atransf = exp,
                      at=log(c(0.05, 0.25, 1, 4, 45)),
                      
                      ilab = cbind(events_short, total_short,
                                   events_long, total_long),
                      
                      ilab.xpos = c(-11.5, -9.5, -6.5, -4.5),
                      
                      xlim = c(-20, 11.17),
                      rows = c(1:7, 13:15),
                      ylim = c(-2, 20),
                      
                      header="Study",
                      pch = 19)
    )
  
  text(c(-11.5, -9.5, -6.5, -4.5),
       x$ylim[2] -1 ,
       rep(c("Events", "Total")),
       cex = 0.9
       )
  
  segments(x0 = -12.5, x1 = -8.7,
           y0 = x$ylim[2] - 0.5, y1 = x$ylim[2] - 0.5)
  
  segments(x0 = -7.5, x1 = -3.7,
           y0 = x$ylim[2] - 0.5, y1 = x$ylim[2] - 0.5)
  
  text(c((-11.5 - 9.5)/2,
         (-6.5 - 4.5)/2),
       x$ylim[2],
       rep(c("Short-Course", "Long-Course")),
       font = 2,
       cex = 0.9, 
       #pos = 4
  )
  
  text(c(-1,1), # X axis
       x$ylim[2] - 3.2, # Y axis
       cex = 0.7, # size
       c("Favors\nShort-Course","Favors\nLong-Course"),
       font = 1,
       pos=c(2,4), # Right + Left aligned
       offset=-1)
  
  text(x$xlim[1],
       pos = 4, # left-aligned
       c(16.5, 8.5), c("RCT","Observational"),
       font = 2)
  
  abline(h=9.5) # horizontal line
  
  RCT_hdi = subset(mortality_hdi, name == "RCT")
  
  metafor::addpoly(x = RCT_hdi$value,
                   ci.lb = RCT_hdi$.lower,
                   ci.ub = RCT_hdi$.upper,
                   
                   mlab = "Total [95% CrI]", # label
                   row=11, # location
                   atransf=exp,
                   efac = 2, # polygon size
                   col=pal[3], border=pal[3]
  )
  
  obs_hdi = subset(mortality_hdi, name == "obs")
  
  metafor::addpoly(x = obs_hdi$value,
                   ci.lb = obs_hdi$.lower,
                   ci.ub = obs_hdi$.upper,
                   
                   mlab = "Total [95% CrI]", # label
                   row=-1, # location
                   atransf=exp,
                   efac = 2, # polygon size
                   col=pal[7], border=pal[7]
  )
  
  
  
  text(x$xlim[1],
       x$ylim[2] + 0.3,
       pos = 4, # left-aligned
       "Mortality",
       cex = 1.5,
       font = 2)
  
}

forest_relapse = function(){
  
  dd = subset(d_logOR, outcome == "relapse30") |> 
    dplyr::arrange(desc(yi))
  
  x = 
    base::with(
      dd,
      
      metafor::forest(x = yi,
                      sei = sei,
                      slab = study,
                      order = desc(type),
                      atransf = exp,
                      at=log(c(0.05, 0.25, 1, 4, 10, 45)),
                      
                      
                      ilab = cbind(events_short, total_short,
                                   events_long, total_long),
                      
                      ilab.xpos = c(-12.5, -10.5, -7.5, -5.5),
                      
                      xlim = c(-22,  12),
                      rows = c(3:7, 14:15),
                      ylim = c(-2, 20),
                      header="Study",
                      pch = 19)
    )
  
  text(c(-12.5, -10.5, -7.5, -5.5),
       x$ylim[2] -1 ,
       rep(c("Events", "Total")),
       cex = 0.9
  )
  
  segments(x0 = -13.5, x1 = -9.7,
           y0 = x$ylim[2] - 0.5, y1 = x$ylim[2] - 0.5)
  
  segments(x0 = -8.5, x1 = -4.7,
           y0 = x$ylim[2] - 0.5, y1 = x$ylim[2] - 0.5)
  
  text(c((-12.5 - 10.5)/2,
         (-7.5 - 5.5)/2),
       x$ylim[2],
       rep(c("Short-Course", "Long-Course")),
       font = 2,
       cex = 0.9, 
       #pos = 4
  )
  
  text(c(-1.25,1.25), # X axis
       x$ylim[2] - 3.2, # Y axis
       cex = 0.7, # size
       c("Favors\nShort-Course","Favors\nLong-Course"),
       font = 1,
       pos=c(2,4), # Right + Left aligned
       offset=-1)
  
  text(x$xlim[1],
       pos = 4, # left-aligned
       c(16.5, 8.5), c("RCT","Observational"),
       font = 2)
  
  abline(h=9.5) # horizontal line
  
  RCT_hdi = subset(relapse_hdi, name == "RCT")
  
  metafor::addpoly(x = RCT_hdi$value,
                   ci.lb = RCT_hdi$.lower,
                   ci.ub = RCT_hdi$.upper,
                   
                   mlab = "Total [95% CrI]", # label
                   row=11, # location
                   atransf=exp,
                   efac = 2, # polygon size
                   col=pal[3], border=pal[3]
  )
  
  obs_hdi = subset(relapse_hdi, name == "obs")
  
  metafor::addpoly(x = obs_hdi$value,
                   ci.lb = obs_hdi$.lower,
                   ci.ub = obs_hdi$.upper,
                   
                   mlab = "Total [95% CrI]", # label
                   row=-1, # location
                   atransf=exp,
                   efac = 2, # polygon size
                   col=pal[7], border=pal[7]
  )
  
  text(x$xlim[1],
       x$ylim[2] + 0.3,
       pos = 4, # left-aligned
       "Relapse",
       cex = 1.5,
       font = 2)
  
  
}


pdf(file = here::here("figures", 
                      'forest_plots.pdf'),
    height=6.5, width=16)

par(mfrow = c(1, 2), # 2 figures, side by side ("1 row, 2 columns")
    # Margins
    mar=c(6, # bottom
          4, # left
          3, # top
          3) # right
) 

forest_mortality()
forest_relapse()

dev.off()


# Density plots ----

theme_set(
  tidybayes::theme_ggdist() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.ticks.y = element_blank(),
      legend.position = 'none',
      plot.margin = margin(20, 20, 20, 20)
    ) 
)

## Mortality ----

prob_mortality_OR = 
  mortality_draws |> 
  dplyr::filter(name %in% c("RCT", "obs")) |> 
  dplyr::group_by(name) |> 
  dplyr::summarise(prob = paste0(round(100*mean(value > 0),1), "%"))

p1 = 
  mortality_draws |> 
  dplyr::filter(name %in% c("RCT", "obs")) |> 
  dplyr::mutate(name = ifelse(name == "RCT", "RCT", "Observational")) |> 
  ggplot() +
  aes(x = value, y = name, fill = name, fill_ramp = stat(x > 0)) +
  tidybayes::stat_halfeye(.width = 0.95,
                          point_interval = tidybayes::median_hdi) +
  ggdist::scale_fill_ramp_discrete(from = "gray85", range = c(0,1)) +
  scale_fill_manual(values = c(pal[7], pal[3])) +
  
  annotate("text", x = log(2), y = 2.6,
           label = prob_mortality_OR[[2,2]], size = 5) +
  
  annotate("text", x = log(2), y = 1.6,
           label = prob_mortality_OR[[1,2]], size = 5) +
  
  
  
  geom_vline(xintercept = 0, linetype = 2, size = 0.6) +
  labs(
    x = "Odds Ratio",
    y = NULL,
    title = "\nMortality"
  ) +
  scale_x_continuous(breaks = log(seq(from = 0.5, to = 2, 0.5)),
                     labels = seq(from = 0.5, to = 2, 0.5)) +
  coord_cartesian(xlim = log(c(0.35, 2.5))) +
  scale_y_discrete(expand = c(0, 0.1))


## Relapse ----

prob_relapse_OR = 
  relapse_draws |> 
  dplyr::filter(name %in% c("RCT", "obs")) |> 
  dplyr::group_by(name) |> 
  dplyr::summarise(prob = paste0(round(100*mean(value > 0),1), "%"))

p2 = 
  relapse_draws |> 
  dplyr::filter(name %in% c("RCT", "obs")) |> 
  dplyr::mutate(name = ifelse(name == "RCT", "RCT", "Observational")) |> 
  ggplot() +
  aes(x = value, y = name, fill = name, fill_ramp = stat(x > 0)) +
  tidybayes::stat_halfeye(.width = 0.95,
                          point_interval = tidybayes::median_hdi) +
  ggdist::scale_fill_ramp_discrete(from = "gray85", range = c(0,1)) +
  scale_fill_manual(values = c(pal[7], pal[3])) +
  
  annotate("text", x = log(3), y = 2.6,
           label = prob_relapse_OR[[2,2]], size = 5) +
  
  annotate("text", x = log(3), y = 1.6,
           label = prob_relapse_OR[[1,2]], size = 5) +
  
  
  
  geom_vline(xintercept = 0, linetype = 2, size = 0.6) +
  labs(
    x = "Odds Ratio",
    y = NULL,
    title = "\nRelapse"
  ) +
  scale_x_continuous(breaks = log(c(seq(from = 0.5, to = 2, 0.5), 3)),
                     labels = c(seq(from = 0.5, to = 2, 0.5), 3)) +
  coord_cartesian(xlim = log(c(0.25, 4))) +
  scale_y_discrete(expand = c(0, 0.1))



## Both ----

(p1 + p2) + patchwork::plot_annotation(tag_levels = "A")

ggsave(width = 10,
       height = 4,
       dpi = 300,
       here("figures", # File path
            "posterior_densities.pdf")) # File name
