# Global constants and functions for
# "Weekly deaths predictive accuracy"
# Jonas Schöley

# Global constants ------------------------------------------------

library(ggplot2)

glob <- list()

glob <- within(glob, {

  # definition of seasons in iso-weeks
  seasons <-
    list(
      northern = list(
        `Winter` = c(49:53, 1:9),
        `Spring` = 10:22,
        `Summer` = 23:35,
        `Fall` = 36:48
      ),
      southern =
        list(
          `Winter` = 23:35,
          `Spring` = 36:48,
          `Summer` = c(49:53, 1:9),
          `Fall` = 10:22
        )
    )
  # iso-week which starts epi-year
  week_epi_year_starts <- 27
  
  # iso-weeks which define usual flu season
  flu_seasons <-
    list(
      northern = c(49:53, 1:12),
      southern = 23:38
    )
  
  # codebook
  codebook <- list(
    sex = c(`m` = "Male", `f` = "Female", `b` = "Total"),
    age_group =
      c(
        `0_14` = "[0,15)",
        `15_64` = "[15,65)",
        `65_74` = "[65,75)",
        `75_84` = "[75,85)",
        `85p` = "85+"
      )
  )
  
  # color coding
  colors <- list(
    sample =
      c(
        training = "grey70",
        test = "grey30"
      ),
    sex =
      c(
        `Male` = "#c60c30",
        `Female` = "#004B87"
      ),
    # color palette
    discrete =
      c('#D23737', # red
        '#3191C9', # blue
        '#D2BC2D', # yellow
        '#4EC93B', # green
        '#881F93', # purple
        '#C5752B') ,# orange
    discrete_light =
      c('#FCB3B3', # red
        '#A7DDFC', # blue
        '#FAEC8E'  # yellow
      )
  )

  MyGGplotTheme <-
    function (
      size = 8,
      family = 'sans',
      scaler = 1,
      axis = 'x',
      panel_border = FALSE,
      grid = 'y',
      show_legend = TRUE,
      ar = NA
    ) {
      
      size_med = size*scaler
      size_sml = round(size*0.7)*scaler
      base_linesize = 0.3*scaler
      
      list(
        theme_minimal(base_size = size_med, base_family = family),
        theme(
          # basic
          text = element_text(color = 'black'),
          line = element_line(size = base_linesize, lineend = 'square'),
          # axis
          #axis.line.y = element_blank(),
          axis.title = element_text(size = size_med, face = 'bold'),
          #axis.ticks = element_line(size = rel(0.5), color = 'black'),
          axis.text = element_text(size = size_med, color = 'black'),
          # strips
          strip.text = element_text(color = 'black', size = size_med),
          strip.background = element_blank(),
          # plot
          title = element_text(face = 'bold'),
          plot.subtitle = element_text(color = 'black', size = size_med, face = 'bold'),
          plot.caption = element_text(color = 'black', size = size_sml, face = 'plain'),
          plot.background = element_blank(),
          panel.background = element_blank(),
          #plot.margin = unit(c(1, 0.1, 0.5, 0.5), units = 'mm'),
          # grid
          panel.grid = element_blank()
        ),
        if (identical(grid, 'y')) {
          theme(panel.grid.major.y =
                  element_line(linetype = 3, color = 'grey80'))
        },
        if (identical(grid, 'x')) {
          theme(panel.grid.major.x =
                  element_line(linetype = 3, color = 'grey80'))
        },
        if (identical(grid, 'xy') | identical(grid, 'yx')) {
          theme(panel.grid.major.y =
                  element_line(linetype = 3, color = 'grey80'),
                panel.grid.major.x =
                  element_line(linetype = 3, color = 'grey80'))
        },
        if (isTRUE(panel_border)) {
          theme(
            panel.border =
              element_rect(fill = NA)
          )
        },
        if (!isTRUE(show_legend)) {
          theme(legend.position = 'none')
        },
        if (axis == 'x') {
          theme(
            axis.line.x = element_line(linetype = 1, color = 'black')
          )
        },
        if (axis == 'y') {
          theme(
            axis.line.y = element_line(linetype = 1, color = 'black')
          )
        },
        if (axis == 'xy') {
          theme(
            axis.line = element_line(linetype = 1, color = 'black')
          )
        },
        if (!is.na(ar)) {
          theme(
            aspect.ratio = ar
          )
        }
      )
    }
  
  fig_dims <- list(
    # figure width (mm)
    width = 170
  )
  
})

# Global functions ------------------------------------------------

#' Convert Week of Year to Date
#'
#' @param year Year integer.
#' @param week Week of year integer (1 to 53).
#' @param weekday Weekday integer (1, Monday to 7, Sunday).
#' @param offset Integer offset added to `week` before date calculation.
#'
#' @return A date object.
#'
#' @source https://en.wikipedia.org/wiki/ISO_8601
#'
#' @author Jonas Schöley
#'
#' @examples
#' # the first Week of 2020 actually starts Monday, December 30th 2019
#' ISOWeekDate2Date(2020, 1, 1)
ISOWeekDate2Date <- function(year, week, weekday = 1, offset = 0) {
  require(ISOweek)
  isoweek_string <-
    paste0(
      year, "-W",
      formatC(
        week + offset,
        flag = "0",
        format = "d",
        digits = 1
      ),
      "-", weekday
    )
  ISOweek2date(isoweek_string)
}

#' Calculate Weeks Since Some Origin Date
#'
#' @param date Date string.
#' @param origin_date Date string.
#' @param week_format Either 'integer' for completed weeks or
#' 'fractional' for completed fractional weeks.
#'
#' @return Time difference in weeks.
#'
#' @author Jonas Schöley
#'
#' @examples
#' # My age in completed weeks
#' WeeksSinceOrigin(Sys.Date(), '1987-07-03')
WeeksSinceOrigin <-
  function(date, origin_date, week_format = "integer") {
    require(ISOweek)
    fractional_weeks_since_origin <-
      as.double(difftime(
        as.Date(date),
        as.Date(origin_date),
        units = "weeks"
      ))
    switch(
      week_format,
      fractional = fractional_weeks_since_origin,
      integer = as.integer(fractional_weeks_since_origin)
    )
  }

#' Return a Sequence of Epi-Year Strings
#' 
#' @param from First year of Epi-Year sequence.
#' @param to Last year of Epi-Year sequence.
#' @param what Type of return value. One of 'slash', 'beginning', 'end'.
#' 
#' @return Sequence of Epi-Years.
#' 
#' @author Jonas Schöley
#' 
#' @examples
#' EpiYearSequence(2000, 2005)
#' # beginning of epi-year
#' EpiYearSequence(2000, 2005, 'beginning')
#' # end of epi-year
#' EpiYearSequence(2000, 2005, 'end')
EpiYearSequence <- function(from, to, what = 'slash') {
  years <- from:to
  switch(
    what,
    slash = paste0(head(years, -1), "/", years[-1]),
    beginning = head(years, -1),
    end = years[-1]
  )
}

# This model and estimates the average mortality rate over
# some years within each week and stratum. The associated
# predict() method multiplies this average mortality with
# given exposures to derive death counts.

AverageMortalityModel <-
  function(df, week, deaths, exposures, ...) {
    require(dplyr)
    .strata <- enquos(...)
    .week <- enquo(week)
    .deaths <- enquo(deaths)
    .exposures <- enquo(exposures)
    
    avg_mx <-
      df %>%
      group_by(!!!.strata, !!.week) %>%
      summarise(
        avg_mortality = mean(!!.deaths / !!.exposures),
        .groups = "drop"
      )
    
    structure(list(avg = avg_mx), class = "avgmx")
  }

predict.avgmx <- function(object, newdata, ...) {
  require(dplyr)
  
  suppressMessages(left_join(newdata, object$avg)) %>%
    pull(avg_mortality)
}

ExportFigure <-
  function(figure,
           path,
           filename,
           width = 170,
           height = 100,
           scale = 1,
           device = 'png',
           dpi = 300,
           add_date = FALSE) {
    require(ggplot2)
    
    if (missing(filename)) {
      filename <- tolower(gsub('\\.', '_', make.names(deparse(substitute(figure)))))
    }
    if (isTRUE(add_date)) {
      filename <- paste0(Sys.Date(), '-', filename)
    }
    
    arguments <-
      list(
        filename = paste0(filename, '.', device),
        plot = figure,
        path = path,
        width = width,
        height = height,
        units = "mm",
        scale = scale,
        dpi = dpi,
        device = device
      )
    if (device == 'pdf') {
      arguments$useDingbats <- FALSE 
    }
    
    do.call(ggsave, arguments)
  }

ExportFiguresFromList <- function(lst, path, ...) {
  figure_names <- tolower(gsub('\\.+', '_', make.names(names(lst))))
  Fun <- function (figure, filename, ...) {
    ExportFigure(figure = figure, filename = filename, ...)
  }
  purrr::pwalk(
    list(lst, figure_names),
    Fun, path = path, ...
  )
}
