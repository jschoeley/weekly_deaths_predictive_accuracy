# Global constants and functions for
# "Weekly deaths predictive accuracy"
#
# 2020-10-19
#
# Jonas Schöley

# Global constants ------------------------------------------------

library(ggplot2)

glob <- list()

glob <- within(glob, {

  countries <- dplyr::tribble(
    ~country_code   , ~country_name              , ~hemisphere,
    'AUT'           , 'Austria'                  , 'northern',
    'BEL'           , 'Belgium'                  , 'northern',
    'BGR'           , 'Bulgaria'                 , 'northern',
    'HRV'           , 'Croatia'                  , 'northern',
    'CZE'           , 'Czech Republic'           , 'northern',
    'DNK'           , 'Denmark'                  , 'northern',
    'GBRTENW'       , 'England And Wales'        , 'northern',
    'EST'           , 'Estonia'                  , 'northern',
    'FIN'           , 'Finland'                  , 'northern',
    'FRATNP'        , 'France'                   , 'northern',
    'DEUTNP'        , 'Germany'                  , 'northern',
    'GRC'           , 'Greece'                   , 'northern',
    'HUN'           , 'Hungary'                  , 'northern',
    'ISL'           , 'Iceland'                  , 'northern',
    'ISR'           , 'Israel'                   , 'northern',
    'ITA'           , 'Italy'                    , 'northern',
    'LVA'           , 'Latvia'                   , 'northern',
    'LTU'           , 'Lithuania'                , 'northern',
    'LUX'           , 'Luxembourg'               , 'northern',
    'NLD'           , 'Netherlands'              , 'northern',
    'GBR_NIR'       , 'Northern Ireland'         , 'northern',
    'NOR'           , 'Norway'                   , 'northern',
    'POL'           , 'Poland'                   , 'northern',
    'PRT'           , 'Portugal'                 , 'northern',
    'RUS'           , 'Russia'                   , 'northern',
    'GBR_SCO'       , 'Scotland'                 , 'northern',
    'SVN'           , 'Slovenia'                 , 'northern',
    'SVK'           , 'Slovakia'                 , 'northern',
    'ESP'           , 'Spain'                    , 'northern',
    'CHE'           , 'Switzerland'              , 'northern',
    'SWE'           , 'Sweden'                   , 'northern',
    'USA'           , 'United States of America' , 'northern'
  )
  
  # iso-week which starts epi-year
  week_epi_year_starts <- 27
  
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

  # iso-weeks which define usual flu season
  flu_seasons <-
    list(
      northern = c(49:53, 1:12),
      southern = 23:38
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
ISOWeekDate2Date <- function (year, week, weekday = 1, offset = 0) {
  require(ISOweek)
  isoweek_string <-
    paste0(
      year, '-W',
      formatC(
        week+offset,
        flag = '0',
        format = 'd',
        digits = 1
      ),
      '-', weekday
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

IsoWeekToEpiWeek <- function (iso_week, w_start = 1, w_53 = FALSE) {
  a <- iso_week - w_start
  ifelse(a < 0, 51 + a + 1 + w_53, a)
}

EpiWeekToIsoWeek <- function (epi_week, w_start = 1, w_53 = FALSE) {
  a <- epi_week + w_start
  ifelse(a > (52 + w_53), a - (52 + w_53), a)
}

#' Export ggplot
#' 
#' @author Jonas Schöley
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

#' Export ggplots Stored in List
#' 
#' @author Jonas Schöley
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

#' Convert population estimates to population exposures
#'
#' Converts population estimates measured at discrete points
#' in time into population exposures by interpolating between
#' the data points using a cubic spline and integrating over
#' arbitrary time intervals.
#' 
#' @param df A data frame.
#' @param x Name of time variable.
#' @param P Name of population variable.
#' @param breaks_out Vector of interpolation points.
#' @param scaler Constant factor to change unit of exposures.
#' @param strata `vars()` specification of variables in df to stratify over.
#'
#' @return A data frame stratified by `strata` with population counts
#' `Px` at time `x1` and exposures `Ex` over time interval `[x1,x2)`. 
#'
#' @author Jonas Schöley, José Manuel Aburto
#'
#' @examples
#' df <-
#'   expand.grid(
#'     sex = c('Male', 'Female'),
#'     age = c('[0, 80)', '80+'),
#'     quarter_week = c(1, 14, 27, 40)
#'   )
#' df$P = rnorm(16, 1e3, sd = 100)
#' Population2Exposures(
#'   df, x = quarter_week, P = P,
#'   breaks_out = 1:57, strata = vars(sex, age)
#' )
Population2Exposures <-
  function (df, x, P, breaks_out, scaler = 1, strata = NA) {
    
    require(dplyr)
    require(tidyr)
    require(purrr)
    
    x = enquo(x); P = enquo(P)
    
    # for each stratum in the data return
    #   - interpolation function
    #   - interpolated population sizes
    #   - interpolated and integrated exposures
    group_by(df, !!!strata) %>% nest() %>%
      mutate(
        # limits of time intervals
        x1 = list(head(breaks_out, -1)),
        x2 = list(breaks_out[-1]),
        # cubic spline interpolation function with linear extrapolation
        interpolation_function = map(
          data,
          ~ splinefun(x = pull(., !!x), y = pull(., !!P),
                      method = 'natural')
        ),
        # interpolated population numbers
        Px = map(
          interpolation_function,
          ~ .(x = head(breaks_out, -1))
        ),
        # interpolated and integrated population exposures
        Ex = map(
          interpolation_function,
          # closed form expression for piecewise polynomial
          # integral exists of course but implementation is
          # left for a later point
          ~ {
            fnct <- .
            map2_dbl(
              unlist(x1),
              unlist(x2),
              ~ integrate(fnct, lower = .x, upper = .y)[['value']]*scaler
            )
          }
        )
      ) %>%
      ungroup() %>%
      select(-data, -interpolation_function) %>%
      unnest(c(x1, x2, Px, Ex))
    
  }
