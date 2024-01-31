#p value labelling function
generate_p_label <- function(f, add_p = TRUE){
  scales::pvalue(f,
                 accuracy = 0.05, # Number to round to
                 decimal.mark = ".", # The character to be used to indicate the numeric decimal point
                 add_p = add_p)
}

#labelling function 
label_assess_id <- function(x){
  output <- lapply(
    paste0(
      x$assess_id, ", ", 
      "Kendall stat: ", round(x$statistic, 2), ", ", 
      generate_p_label(x$p.value)
    ),
    HTML)
  return(output)
}

#colour scale range function
scale_range <- function(x){
  scale_lim <- max(abs(range(x)))
  scale_range <- c(-1*scale_lim, scale_lim)
  return(scale_range)
}

#generate the base map function
generate_map <- function(x, lf, legend=FALSE){
  
  #subset to include only valid statistics
  x <- subset(x, !is.na(statistic))
  
  #generate the colour palette function
  pal <- colorNumeric(
    palette = c("#5e3c99","#ffffbf","#e66101"),
    domain = scale_range(x$statistic)
  )
  
  #title function
  tag.map.title <- tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 12px;
  }
"))
  
  #generate title
  title <- tags$div(
    tag.map.title, HTML(
      str_to_title(
        gsub("_", " ", lf)
      )
    )
  )  
  
  #subset to focus on selected lifeform 
  temp <- x %>%
    filter(lifeform==lf)
  
  #data subset to COMP4 layer
  temp_COMP4 <- temp %>%
    filter(type == "COMP4")
  
  #data subset for WFD layer
  temp_WFD <- temp %>%
    filter(type == "WFD")
  
  #data subset for stations layer
  temp_stations <- temp %>%
    filter(type == "Stations")
  
  #determine the layers available
  temp_layers <- unique(temp$type)
  
  #generate leaflet
  temp_map <- leaflet() %>%
    addControl(title, position = "topleft", className="map-title")
  
  # set max zoom for when only one station to preserve geographic context
  if(length(unique(temp$assess_id)) == 1){
    temp_map <- temp_map %>%
      addProviderTiles(providers$CartoDB.PositronNoLabels,
                       options = providerTileOptions(maxZoom = 5))
  } else {
    temp_map <- temp_map %>%
      addProviderTiles(providers$CartoDB.PositronNoLabels)
  }
  
  if("COMP4" %in% temp_layers){
    temp_map <- temp_map %>%
      addMapPane("COMP4", zIndex = 410) %>%
      addPolygons(data=temp_COMP4,
                  group = "COMP4",
                  options = pathOptions(pane = "COMP4"),
                  layerId=~assess_id,
                  color = ~pal(temp_COMP4$statistic),
                  label = label_assess_id(temp_COMP4),
                  weight = 1, smoothFactor = 0.5, opacity = 1.0, fillOpacity = 0.5
      )
  }
  
  if("WFD" %in% temp_layers){
    temp_map <- temp_map %>%
      addMapPane("WFD", zIndex = 420) %>%
      addPolygons(data=temp_WFD,
                  group = "WFD",
                  options = pathOptions(pane = "WFD"),
                  layerId=~assess_id,
                  color = ~pal(temp_WFD$statistic),
                  label = label_assess_id(temp_WFD),
                  weight = 1, smoothFactor = 0.5, opacity = 1.0, fillOpacity = 0.5
      )
  }
  
  if("Stations" %in% temp_layers){
    temp_map <- temp_map %>%
      addMapPane("Stations", zIndex = 430) %>%
      addCircleMarkers(data=temp_stations,
                       group = "Stations",
                       options = pathOptions(pane = "Stations"),
                       layerId=~assess_id,
                       color = ~pal(temp_stations$statistic),
                       label = label_assess_id(temp_stations),
                       weight = 1, opacity = 1.0, fillOpacity = 0.5
      )
  }
  
  temp_map <- temp_map %>%
    addLayersControl(overlayGroups = temp_layers)
  
  if(legend){
    temp_map <- temp_map %>%
      addLegend("bottomright", pal = pal, values = scale_range(x$statistic),
                title = "Kendall statistic",
                opacity = 1)
  } 
  
  return(temp_map)
}

#slider date range function
slider_date_range <- function(x){
  y <- c(
    as.Date(paste(min(as.numeric(x)), 1, 1, sep="-")),
    as.Date(paste(max(as.numeric(x)), 12, 31, sep="-"))
  )
  return(y)
}

#time-series plotting function
plot_ts <- function(x, lf, assessment_range, comparison_range, text_string) {
  
  temp <- x %>%
    mutate(count_month = round(count, 2)) %>%
    mutate(date_month = as.Date(paste(year, month, 15), "%Y %m %d"),
           date_year = as.Date(paste(year, 7, 2), "%Y %m %d")) %>%
    group_by(date_year) %>%
    mutate(count_year = round(mean(count_month, na.rm=T), 2)) %>%
    ungroup()
  
  y_range <- range(temp$count_month)
  a_range <- slider_date_range(assessment_range)
  c_range <- slider_date_range(comparison_range)
  
  units <- x$abundance_type_units[1]
  
  suppressWarnings({
    ggplot() +
      geom_rect(aes(xmin = a_range[1],
                    xmax = a_range[2],
                    ymin = y_range[1],
                    ymax = y_range[2]), fill = '#48B2DE', alpha=0.2) +
      geom_rect(aes(xmin = c_range[1],
                    xmax = c_range[2],
                    ymin = y_range[1],
                    ymax = y_range[2]), fill = '#E57872', alpha=0.2) +
      geom_line(data = temp, aes(x = date_month, y = count_month), colour="blue", linewidth=0.3) +
      geom_point(data = temp, aes(x = date_month, y = count_month, text = paste0("log10(", lf, " ", units, ")", ": ", count_month, "<br>",
                                                                                 "Month: ", format(as.Date(paste(year, month, 15, sep="-")), "%Y-%m"))),
                 shape = 21, stroke=0.1, size=0.5, fill="blue") + 
      geom_smooth(data = temp, aes(x = date_year, y = count_year), formula = y ~ x,
                  linetype="dashed", colour="red", 
                  method = 'lm', se = FALSE) +
      geom_line(data = temp, aes(x = date_year, y = count_year)) +
      geom_point(data = temp, aes(x = date_year, y = count_year, text = paste0("log10(", lf, " ", units, ")", ": ", count_year, "<br>",
                                                                               "Year: ", format(as.Date(paste(year, 07, 2, sep="-")), "%Y"))),
                 shape = 21, stroke=0.2, size=1.5, fill="grey80") +
      scale_y_continuous(name=paste0("log10(", lf, " ", units, ")")) +
      scale_x_date(limits = slider_date_range(format(temp$date_month, "%Y"))) +
      ggtitle(text_string) +
      theme_minimal() +
      theme(axis.title.x=element_blank(),
            plot.title = element_text(hjust = 0.5,
                                      size = 10))
  })
}

#function for selecting subsets of the data for assessment and comparison periods
dataSelect <- function(x, lims, lf){
  
  df_temp_output <- x %>%
    filter(year >= lims[1],
           year <= lims[2]) %>%
    pivot_wider(names_from = lifeform, values_from = count) %>%
    relocate(assess_id) %>%
    arrange(assess_id, year, month) %>%
    mutate(n = as.integer(n)) %>%
    dplyr::select(unique_id, assess_id, year, month, n, abundance_type_units, all_of(lf))
  
  return(df_temp_output)
}

#quality control steps to ensure PI results are reliable, based on assessment period data
qc_assessment <- function(x, lf, min_months = 10, robust_months = 30, ind_years = 3, rep_months = 1, prop_nonzero = 0.2){
  
  if(!is.null(x)){
    output <- x %>%
      dplyr::select(assess_id, year, month, n, all_of(lf)) %>%
      mutate_at(vars(c(lf[1], lf[2])), function(f) 10^f - 1) %>%
      mutate(lf1_perc_zero = sum(round(get(lf[1])) == 0) / nrow(.)) %>%
      mutate(lf2_perc_zero = sum(round(get(lf[2])) == 0) / nrow(.)) %>%
      mutate(lf_count_min_nonzero = ifelse(mean(lf1_perc_zero) <= prop_nonzero &
                                             mean(lf2_perc_zero) <= prop_nonzero,
                                           TRUE,
                                           FALSE)) %>%
      dplyr::select(-c(lf[1], lf[2], lf1_perc_zero, lf2_perc_zero)) %>%
      mutate(n_months_min_sufficient = length(unique(paste(year,month))) >= min_months) %>% # at least 10 individual months represented for minimal PI
      mutate(n_months_robust_sufficient = length(unique(paste(year,month))) >= robust_months) %>% # at least 30 individual months represented for robust PI
      mutate(n_years_sufficient = length(unique(year)) >= ind_years) %>% # at least three years of non-interpolated samples
      group_by(assess_id, month) %>%
      mutate(n_months_rep_sufficient = length(unique(year)) >= rep_months) %>% # each month represented at least twice
      ungroup() %>%
      mutate(n_months_rep_sufficient = all(n_months_rep_sufficient)) %>%
      rename_with(.fn = ~paste(., prop_nonzero, sep=":"), .cols = lf_count_min_nonzero) %>%
      rename_with(.fn = ~paste(., min_months, sep=":"), .cols = n_months_min_sufficient) %>%
      rename_with(.fn = ~paste(., robust_months, sep=":"), .cols = n_months_robust_sufficient) %>%
      rename_with(.fn = ~paste(., ind_years, sep=":"), .cols = n_years_sufficient) %>%
      rename_with(.fn = ~paste(., rep_months, sep=":"), .cols = n_months_rep_sufficient) %>%
      left_join(x)
    
    return(output)
  } else {
    return(NULL)
  }
}

#quality control steps to ensure PI results are reliable, based on comparison period data and no overlap with assessment
qc_comparison <- function(comp, assess, prop_comp = 0.25, n_overlap = 1){
  
  if(!is.null(comp) &
     !is.null(assess)){
    
    temp_comp <- comp %>%
      dplyr::select(assess_id, year, month, n) %>%
      mutate(type = "comp")
    
    temp_assess <- assess %>%
      dplyr::select(assess_id, year, month, n) %>%
      mutate(type = "assess") %>%
      bind_rows(temp_comp) %>%
      mutate(date = as.Date(paste(year, month, 15, sep="-"))) %>%
      pivot_wider(names_from = type, values_from = date)
    
    #test to see if assessment and comparison data overlap
    overlap_not_present <- if(min(temp_assess$assess, na.rm=T) <= max(temp_assess$comp, na.rm=T) &&
                              min(temp_assess$comp, na.rm=T) <= max(temp_assess$assess, na.rm=T)){
      FALSE
    } else {
      TRUE
    } 
    
    #assess the comparison data to ensure it meets the criteria for generating a reliable PI calculation
    output <- temp_comp %>%
      dplyr::select(assess_id, year, month, n) %>%
      mutate(prop_monitored = n()/(12*length(seq(min(year), max(year), 1)))) %>%
      mutate(prop_comp_monitored_min = ifelse(prop_monitored >= prop_comp,
                                              TRUE,FALSE)) %>%
      dplyr::select(-prop_monitored) %>%
      rename_with(.fn = ~paste(., prop_comp, sep=":"), .cols = prop_comp_monitored_min) %>%
      mutate(overlap_not_present = overlap_not_present) %>%
      rename_with(.fn = ~paste(., n_overlap, sep=":"), .cols = overlap_not_present)
    
    return(output)
  } else {
    return(NULL)
  }
}

#setup the quality control output
qc_text_output <- function(x){
  
  if(!is.null(x)){
    
    #quality control steps for output message
    temp_qc <- x %>%
      dplyr::select(assess_id, 
                    contains(":")) %>%
      distinct() %>%
      pivot_longer(!assess_id, names_to="condition", values_to="tf") %>%
      separate_wider_delim(condition, ":", names=c("condition", "criterion"))
    
    return(temp_qc)
  } else {
    return(NULL)
  }
}

# display debugging messages in R (if local) 
# and in the console log (if running in shiny)
debug_msg <- function(...) {
  is_local <- Sys.getenv('SHINY_PORT') == ""
  in_shiny <- !is.null(shiny::getDefaultReactiveDomain())
  txt <- toString(list(...))
  if (is_local) message(txt)
  if (in_shiny) shinyjs::runjs(sprintf("console.debug(\"%s\")", txt))
}

#function for creating labeller lookup table
generate_ts_label <- function(x, assess_id_label, units_label){
  
  #create labeller lookup table
  plot_label <- paste0(
    assess_id_label, ', ',
    'Kendall stat: ', round(x$statistic, 3), ', ',
    'p: ', generate_p_label(x$p.value, add_p = FALSE), '\n',
    'n: ', x$n, ', ',
    "Sen's slope: ", 
    ifelse(!is.na(x$sens_estimate),
           paste0(round(x$sens_estimate), " ", units_label, "/", "year", ', '),
           paste0(NA)
    ),
    ifelse(!is.na(x$sens_p.value),
           paste0('p: ', 
                  generate_p_label(x$sens_p.value, add_p = FALSE)),
           paste0("")
    )
  )
  return(plot_label)
}

#function for generating the two ts_plots
generate_ts_plot <- function(x, y, lf_select, assessment_range, comparison_range){
  
  df_subset <- x %>%
    filter(lifeform == lf_select)
  
  #generate text string for plot labeling
  temp_ken <- y %>%
    filter(lifeform==lf_select,
           assess_id==df_subset$assess_id[1])
  
  #generate plot label using custom labeller function
  plot_label <- generate_ts_label(temp_ken, 
                                  assess_id_label=df_subset$assess_id[1],
                                  units_label=df_subset$abundance_type_units[1]
  )
  
  #generate ts_plot using custom function
  ts_plot <- plot_ts(df_subset,
                     lf = lf_select,
                     assessment_range = assessment_range,
                     comparison_range = comparison_range,
                     text_string = plot_label)
  
  # Convert ggplot to plotly with ggplotly
  ts_plot_plotly <- ggplotly(
    ts_plot,
    tooltip = "text"  # Use the "text" aesthetic for tooltips
  )
  
  return(ts_plot_plotly)
}

#function for creating plot label
generate_pi_label <- function(x, y, assess_id_label){
  
  #create labeller lookup table
  plot_label <- paste0(assess_id_label, ', ',
                       'PI: ', round(x$PI, 3),'\n',
                       'assess points: ', nrow(y), ',', ' ',
                       'comp points: ', x$newPoints, ',', ' ',
                       'binom-p: ', generate_p_label(x$binomial_probability, add_p = FALSE), ',', '  ',
                       'chi-sq: ', round(x$chi.sq, 3)
  )
  
  return(plot_label)
}

