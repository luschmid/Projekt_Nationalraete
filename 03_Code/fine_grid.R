

theme_bw_finegrid <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.key        = element_rect(colour = "grey80"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major.x =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.minor.x =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.spacing =      unit(0.1, "lines"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50")
    )
}


theme_bw_finegrid_finer <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.key        = element_rect(colour = "grey80"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major.x =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.minor.x =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.spacing =      unit(0.1, "lines"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50")
    )
}


theme_bw_finegrid_special <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.key        = element_rect(colour = "grey80"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major.x =  element_line(colour = "grey", size = 0.5,linetype="solid"),
      panel.grid.minor.x =  element_line(colour = "grey", size = 0.5,linetype="solid"),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.spacing =  element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.margin =      unit(0.1, "lines"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50")
    )
  
}

theme_bw_finegrid_horizontal <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.key        = element_rect(colour = "grey80"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5,linetype="solid"),
      panel.grid.minor.y =   element_blank(),
      panel.spacing =      unit(0.1, "lines"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50")
    )
}



theme_bw_finegrid_horizontal_minor <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.key        = element_rect(colour = "grey80"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5,linetype="solid"),
      panel.grid.minor.y =    element_line(colour = "grey", size = 0.5,linetype="solid"),
      panel.margin =      unit(0.1, "lines"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50")
    )
}

theme_bw_finegrid_horizontal_fine_x <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.key        = element_rect(colour = "grey80"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major.x = element_line(colour = "grey", size = 0.5,linetype="dotted"),
      panel.grid.minor.x =  element_blank(),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5,linetype="solid"),
      panel.grid.minor.y =   element_blank(),
      panel.margin =      unit(0.1, "lines"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50")
    )
}
