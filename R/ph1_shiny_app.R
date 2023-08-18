#install and/or load the required packages
list.of.packages <- c("data.table",
                      "fst",
                      "leaflet",
                      "dplyr",
                      "shiny",
                      "sf",
                      "tidyverse",
                      "rgdal",
                      "broom",
                      "dichromat",
                      "pracma",
                      "shinyjs",
                      "plotly",
                      "shinyWidgets",
                      "trend",
                      "scales",
                      "EnvStats",
                      "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

#enter the filename of the multipolygon shapefile the data is partitioned with
file_shp_part <- "COMP4_WFD.shp"

#enter the filename of the point shapefile the station data is partitioned with
file_pts_part <- "Stations.shp"

#enter the main directory to use to access the processed data
dir_raw <- "../data/"

#generate a string for the directory of the multipolygon shapefile
dir_shp <- paste0(dir_raw, gsub(".shp","",file_shp_part), "/")

#generate a string for the directory of the point shapefile
dir_pts <- paste0(dir_raw, gsub(".shp","",file_pts_part), "/")

#load the supporting functions for calculating the indicator
source("../R/supporting_functions/PI_functions_v1.R")

#read the plankton data
df_plot <- read_fst(path=paste0(dir_raw, "COMP4_WFD_Stations_plankton", ".fst")) %>%
  filter(tconf >= 0.5) %>%
  group_by(data_id, assess_id, lifeform) %>%
  filter(length(unique(year)) >= 3) %>%
  ungroup()
  
#load the shapefile associated with the data
#shp_part <- st_as_sf(rgdal::readOGR(paste(dir_shp, file_shp_part, sep="")))
shp_part <- st_read(paste(dir_shp, file_shp_part, sep=""))

#generate points shape file to represent the stations
#shp_pts <- st_as_sf(rgdal::readOGR(paste(dir_pts, file_pts_part, sep="")))
shp_pts <- st_read(paste(dir_pts, file_pts_part, sep=""))

#merge the points and polygons into a single shapefile
shp_merged <- rbind(shp_part, shp_pts)

#acceptable lifeform pair combinations
#construct a dataframe of relevant lifeform pair comparisons
df_lf <- rbind(data.frame(V1 = "diatom", V2 = "dinoflagellate"),
               data.frame(V1 = "tychopelagic_diatoms", V2 = "pelagic_diatoms"),
               data.frame(V1 = "lg_copepods", V2 = "sm_copepods"),
               data.frame(V1 = "holoplankton", V2 = "meroplankton"),
               data.frame(V1 = "lg_phyto", V2 = "sm_phyto"),
               data.frame(V1 = "phytoplankton", V2 = "non_carniv_zoo"),
               data.frame(V1 = "crustacean", V2 = "gelatinous"),
               data.frame(V1 = "gelatinous", V2 = "fish_larvae")
)

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
plot_ts <- function(x, lf, text_string) {
  
  temp <- x %>%
    mutate(count_month = round(count, 2)) %>%
    mutate(date_month = as.Date(paste(year, month, 15), "%Y %m %d"),
           date_year = as.Date(paste(year, 7, 2), "%Y %m %d")) %>%
    group_by(date_year) %>%
    mutate(count_year = round(mean(count_month, na.rm=T), 2)) %>%
    ungroup()
  
  ggplot() +
    geom_line(data = temp, aes(x = date_month, y = count_month), colour="blue", linewidth=0.3) +
    geom_point(data = temp, aes(x = date_month, y = count_month),
               shape = 21, stroke=0.1, size=0.5, fill="grey80") + 
    geom_smooth(data = temp, aes(x = date_year, y = count_year), formula = y ~ x,
                linetype="dashed", colour="red", 
                method = 'lm', se = FALSE) +
    geom_line(data = temp, aes(x = date_year, y = count_year)) +
    geom_point(data = temp, aes(x = date_year, y = count_year), shape = 21, stroke=0.2, size=1.5, fill="grey80") +
    scale_y_continuous(name=paste0("log10(", lf, " ", temp$abundance_type_units[1], ")")) +
    scale_x_date(limits = c(min(temp$date_month), max(temp$date_month))) +
    ggtitle(text_string) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    size = 10))
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

#quality control steps to ensure PI results are reliable, based on comparison period data
qc_comparison <- function(x, prop_comp = 0.25){
  
  if(!is.null(x)){
  
    #assess the comparison data to ensure it meets the criteria for generating a reliable PI calculation
    output <- x %>%
      dplyr::select(assess_id, year, month, n) %>%
      mutate(prop_monitored = n()/(12*length(seq(min(year), max(year), 1)))) %>%
      mutate(prop_comp_monitored_min = ifelse(prop_monitored >= prop_comp,
                                              TRUE,FALSE)) %>%
      rename_with(.fn = ~paste(., prop_comp, sep=":"), .cols = prop_comp_monitored_min) %>%
      dplyr::select(-prop_monitored) %>%
      left_join(x)
    
    return(output)
  } else {
    return(NULL)
  }
}

#setup the quality control output
qc_output <- function(x){
  
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

##################################################################################################################################
# UI
##################################################################################################################################
ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs features
  
  titlePanel("PH1/FW5 Explorer"),
  
  # Add a text output to display test data
  #tableOutput("test"),

  sidebarLayout(
    div(
      # Wrap sidebarPanel in a div with class "sidebar-panel"
      sidebarPanel(
        selectInput('dataset', 'Select dataset', choices = c("Please make a selection", unique(df_plot$data_name))),
        selectInput('lifeform_pair', 'Select lifeform pair', choices = NULL),  # Empty initially
        selectInput('assessment_unit', 'Select assessment unit or click one on map', choices = NULL),  # Empty initially
        # Input: Specification of range within an interval ----
        setSliderColor(c("blue", "#48B2DE", "#E57872"), c(1, 2, 3)),
        sliderInput("range_slider", "Year range for calculating time-series trend:",
                    min = 1900, max = as.numeric(format(Sys.Date(), "%Y")),
                    value = c(1900, as.numeric(format(Sys.Date(), "%Y"))),
                    sep = "",
                    step = 1,
                    dragRange = TRUE),
        # New sliders for Assessment and Comparison
        sliderInput("assessment_slider", "Plankton Index assessment period (envelope):", 
                    min = 1900, max = as.numeric(format(Sys.Date(), "%Y")),
                    value = c(1900, as.numeric(format(Sys.Date(), "%Y"))),
                    sep = "",
                    step = 1,
                    dragRange = TRUE),
        sliderInput("comparison_slider", "Plankton Index comparison period (points):", 
                    min = 1900, max = as.numeric(format(Sys.Date(), "%Y")),
                    value = c(1900, as.numeric(format(Sys.Date(), "%Y"))),
                    sep = "",
                    step = 1,
                    dragRange = TRUE),
        # Add a placeholder div for the error message
        tags$div(id = "error_message", style = "color: red; font-weight: bold;"),
        
        # Conditional panel to show env_plot when a dataset is selected
        conditionalPanel(
          condition = "input.dataset != null",
          id = "env_plot_conditional",  # Add id attribute
          plotlyOutput("env_plot", height = "500px")
        ),
        
        # Conditional panel to show env_text when df_pi() is a character string
        conditionalPanel(
          condition = "typeof Shiny.outputBindings.bindingNames['shiny.textOutput'] !== 'undefined' && 
                        Shiny.outputBindings.bindingNames['shiny.textOutput'].includes('env_text')",
          textOutput("env_text")
        )
      ),
      class = "sidebar-panel" # Add the "sidebar-panel" class to the div
    ),
    # Main panel with leafletOutput and ts_plot
    mainPanel(
      fluidRow(
        column(width = 6, leafletOutput("map1")),
        column(width = 6, plotlyOutput("ts_plot1"), style = "height: 100%")  # Set the height of plotlyOutput to 100%
      ),
      fluidRow(
        column(width = 6, leafletOutput("map2")),
        column(width = 6, plotlyOutput("ts_plot2"), style = "height: 100%")  # Set the height of plotlyOutput to 100%
      ),
      fluidRow(
        column(width = 6, htmlOutput("user_prompt")),  # Empty column to align left edge
        column(width = 6, plotlyOutput("histo_plot", height = "200px"))  # Adjust maxHeight as needed
      )
    )
  )
)

##################################################################################################################################
# SERVER
##################################################################################################################################

server <- function(input, output, session) {

  # Create a reactive value to track whether a dataset has been selected
  dataset_selected <- reactiveVal(FALSE)

  # Create separate reactive values for storing the clicked IDs for map1 and map2
  clicked_id <- reactiveVal(NULL)
  
  # Create a reactive value to track whether a lifeform pair is selected
  lifeform_selected <- reactiveVal(FALSE)
  
  #create empty reactive vals for the filtered dataframes
  df_dset <- reactiveVal(NULL)
  df_dset_range <- reactiveVal(NULL)
  df_dset_range_lf <- reactiveVal(NULL)
  df_dset_range_lf_aid <- reactiveVal(NULL)
  
  # Function to update the error message
  updateErrorMessage <- function(message) {
    shinyjs::runjs(sprintf("$('#error_message').text('%s');", message))
  }
  
  # Observe event to update error message for ts_plot1
  observeEvent(input$ts_plot1, {
    debug_msg("ts_plot1 NULL")
    if (is.null(input$ts_plot1)) {
      updateErrorMessage("Please select a map unit to generate plot output")
    } else {
      updateErrorMessage("")
    }
  })
  
  # Observe event to update error message for ts_plot2
  observeEvent(input$ts_plot2, {
    debug_msg("ts_plot2 NULL")
    if (is.null(input$ts_plot2)) {
      updateErrorMessage("Please select a map unit to generate plot output")
    } else {
      updateErrorMessage("")
    }
  })
  
  # Observe event to update error message for env_plot
  observeEvent(input$env_plot, {
    debug_msg("env_plot NULL")
    if (is.null(input$env_plot)) {
      updateErrorMessage("Please select a map unit to generate plot output")
    } else {
      updateErrorMessage("")
    }
  })
  
  # Update the lifeform_pair and assessment_unit dropdowns when the dataset is selected
  observeEvent(input$dataset, {
    req(input$dataset)
    debug_msg("Waiting for user to select a dataset before lifeform and assessment area options generated")
    
    #reset clicked mapid to NULL
    clicked_id(NULL)
    
    # Reset the reactive objects for filtered_data
    lifeform_selected(FALSE)
    
    # Construct valid_choices for lifeform_pair dropdown
    lifeform_pair_choices <- do.call(paste, c(df_lf[df_lf$V1 %in% unique(df_plot$lifeform[df_plot$data_name == input$dataset]) &
                                              df_lf$V2 %in% unique(df_plot$lifeform[df_plot$data_name == input$dataset]),], sep = "-"))
    
    # Construct valid_choice_labels for assessment_unit dropdown
    assessment_unit_choices <- sort(unique(df_plot$assess_id[df_plot$data_name == input$dataset]))
    
    # Update the lifeform_pair and assessment_unit dropdowns
    updateSelectInput(session, 'lifeform_pair', choices = c("Please make a selection", lifeform_pair_choices))
    updateSelectInput(session, 'assessment_unit', choices = c("Please make a selection", assessment_unit_choices))
    
    # Update the dataset_selected reactive value to indicate that a dataset has been selected
    if(input$dataset != "Please make a selection"){
      dataset_selected(TRUE)
    } else {
      dataset_selected(FALSE)
    }
    
    # Check if there's only one option available for lifeform_pair and assessment_unit
    if (length(lifeform_pair_choices) == 1) {
      updateSelectInput(session, 'lifeform_pair', selected = lifeform_pair_choices)
    }
    
    if (length(assessment_unit_choices) == 1) {
      updateSelectInput(session, 'assessment_unit', selected = assessment_unit_choices)
      clicked_id(assessment_unit_choices)
    }
  })
  

  # Observe the dataset_selected reactive value to update the UI elements based on its value
  observe({
    if (dataset_selected()) {
      shinyjs::show("lifeform_pair")
      shinyjs::show("assessment_unit")
    } else {
      shinyjs::hide("lifeform_pair")
      shinyjs::hide("assessment_unit")
    }
  })
  
  # Observe the lifeform_selected reactive value to update the UI elements based on its value
  observe({
    if (lifeform_selected()) {
      shinyjs::show("map1")
      shinyjs::show("map2")
    } else {
      shinyjs::hide("map1")
      shinyjs::hide("map2")
    }
  })
  
  # Update the lifeform_selected reactive value when a lifeform pair is selected
  observe({
    if (!is.null(input$lifeform_pair) &
        input$lifeform_pair != "Please make a selection") {
      lifeform_selected(TRUE)
    } else {
      lifeform_selected(FALSE)
    }
  })
  
  # Split the lifeform pairs string into the two lifeforms
  lifeforms <- reactive({
    if(!is.null(input$lifeform_pair) &
      input$lifeform_pair != "Please make a selection"){
      output <- str_split(input$lifeform_pair, "-")[[1]]

      return(output)
    } else {
      return(NULL)
    }
  })
  
  #filter the main dataframe by dataset
  df_dset <- reactive({
    if(!is.null(input$dataset) &
       input$dataset != "Please make a selection"){
      
      debug_msg("Filtering on dataset input")
      
      # Filter the data by dataset and lifeform pair
      df_temp_output <- df_plot %>%
        filter(data_name == input$dataset)
      
      if(nrow(df_temp_output) > 0){
        
        df_temp_output <- df_temp_output %>%
          mutate(year = as.integer(year),
                 month = as.integer(month),
                 date = as.Date(paste(year, month, 15, sep = "-")))
        
      } else {
        df_temp_output <- NULL
      }
    } else {
      df_temp_output <- NULL
    }
    return(df_temp_output)
  })
  
  # Update slider input depending on the selected dataset
  observeEvent(input$dataset, {
    req(df_dset())
    debug_msg("Range slider input updated")
    
    updateSliderInput(session, 'range_slider',
                      min = min(df_dset()$year),
                      max = max(df_dset()$year),
                      value = c(min(df_dset()$year), max(df_dset()$year))
    )
  })
  
  #filter the dataframe that has already been filtered by dataset and lifeform pair by range_slider
  df_dset_range <- reactive({
    if(!is.null(df_dset()) &
       !is.null(input$range_slider)){
      debug_msg("Filtering data by range slider")
      
      df_temp_output <- df_dset() %>%
        filter(year >= min(input$range_slider),
               year <= max(input$range_slider))
    } else {
      df_temp_output <- NULL
    }
    return(df_temp_output)
  })

  #filter the dataframe that has already been filtered by dataset by lifeform pair
  df_dset_range_lf <- reactive({
    if(!is.null(df_dset_range()) &
       !is.null(lifeforms())){
      debug_msg("Filtering data by lifeform pair")
      
      df_temp_output <- df_dset_range() %>%
        filter(lifeform == lifeforms()[1] | lifeform == lifeforms()[2])
    } else {
      df_temp_output <- NULL
    }
    return(df_temp_output)
  })
  
  
  
  #Kendall trend test and Sen's slope calculation
  df_ken <- reactive({
    
    if(!is.null(df_dset_range_lf())){
      debug_msg("Kendall trend test and Sen's slope calculations")
      
      df_temp_output <- df_dset_range_lf() %>%
        group_by(lifeform, assess_id, is_point, year) %>%
        summarise(count_year = mean(count, na.rm = TRUE),
                  n = sum(n, na.rm=T),
                  .groups='drop') %>%
        group_by(lifeform, assess_id) %>%
        filter(n() >= 3) %>%
        ungroup()
      
      if(nrow(df_temp_output) > 0){
        
        df_temp_output <- df_temp_output %>%
          mutate(year = as.integer(year)) %>%
          group_by(lifeform, assess_id, is_point) %>%
          nest() %>%
          mutate(n = map_dbl(data, ~sum(.x$n, na.rm=T))) %>%
          mutate(fits = map(data, ~EnvStats::kendallTrendTest(count_year ~ year, ci.slope = FALSE, data = .x)),
                 fits2 = map(fits, ~structure(.x, class = "htest")),
                 fits3 = map(fits2, ~tidy(.x))) %>%
          dplyr::select(-c(fits, fits2)) %>%
          unnest(cols=c(fits3)) %>%
          mutate(data = map(data, ~mutate(.x, count_year_lin = 10^count_year - 1))) %>% # regenerate linear count to calculate Sen's slope
          mutate(sens_estimate = ifelse(p.value <= 0.05, 
                                        map_dbl(data, ~trend::sens.slope(.x$count_year_lin)$estimates), 
                                        NA),
                 sens_p.value = ifelse(p.value <= 0.05, 
                                       map_dbl(data, ~trend::sens.slope(.x$count_year)$estimates), 
                                       NA)) %>%
          mutate(sens_estimate = as.numeric(sens_estimate),
                 sens_p.value = as.numeric(sens_p.value)) %>%
          dplyr::select(lifeform, assess_id, is_point, n, statistic, p.value, sens_estimate, sens_p.value)
        
      } else {
        
        df_temp_output <- NULL
      }
      
    } else {
      
      df_temp_output <- NULL
      
    }
    
    return(df_temp_output)
    
  })
  
  # Generate map1
  output$map1 <- renderLeaflet({
    if (lifeform_selected() &
        !is.null(df_ken()) &
        !is.null(lifeforms())) {
      debug_msg("Generating map1")
      
      output <- generate_map(merge(shp_merged, df_ken(), by = c("assess_id", "is_point")), lf = lifeforms()[1])
    } else {
      output <- NULL
    }
    return(output)
  })
  
  # Generate map2
  output$map2 <- renderLeaflet({
    if (lifeform_selected() &
        !is.null(df_ken()) &
        !is.null(lifeforms())) {
      debug_msg("Generating map2")
      
      output <- generate_map(merge(shp_merged, df_ken(), by = c("assess_id", "is_point")), lf = lifeforms()[2], legend = TRUE)
    } else {
      output <- NULL
    }
    return(output)
  })
  
  # Create reactive values to store the current bounds of both maps
  map1_bounds <- reactiveVal()
  map2_bounds <- reactiveVal()
  
  # Add a delay flag and a timer to reset it
  delay_update <- reactiveVal(FALSE)
  delay_timer <- reactiveTimer(200)  # 200 milliseconds
  
  # Function to update map1 bounds (called when map1 bounds change)
  updateMap1Bounds <- function(mapBounds) {
    if (!delay_update()) {
      map1_bounds(mapBounds)
      delay_update(TRUE)
      delay_timer()  # Start the timer
      # Update map2 only if the bounds change is triggered by the user, not automatic fitting
      if (!identical(mapBounds, map2_bounds())) {
        map2_bounds(mapBounds)
        leafletProxy("map2") %>%
          fitBounds(mapBounds$west, mapBounds$south, mapBounds$east, mapBounds$north)
      }
    }
  }
  
  # Function to update map2 bounds (called when map2 bounds change)
  updateMap2Bounds <- function(mapBounds) {
    if (!delay_update()) {
      map2_bounds(mapBounds)
      delay_update(TRUE)
      delay_timer()  # Start the timer
      # Update map1 only if the bounds change is triggered by the user, not automatic fitting
      if (!identical(mapBounds, map1_bounds())) {
        map1_bounds(mapBounds)
        leafletProxy("map1") %>%
          fitBounds(mapBounds$west, mapBounds$south, mapBounds$east, mapBounds$north)
      }
    }
  }
  
  # Observe event when the delay timer expires and reset the delay flag
  observeEvent(delay_timer(), {
    delay_update(FALSE)
  })
  
  # Observe event when map1 bounds change and update map2
  observeEvent(input$map1_bounds, {
    debug_msg("Updating map1 bounds")
    updateMap1Bounds(input$map1_bounds)
  })
  
  # Observe event when map2 bounds change and update map1
  observeEvent(input$map2_bounds, {
    debug_msg("Updating map2 bounds")
    updateMap2Bounds(input$map2_bounds)
  })
  
  # Function to update the clicked ID
  observeEvent(input$map1_shape_click, {
    clicked_id(input$map1_shape_click$id)
  })
  observeEvent(input$map2_shape_click, {
    clicked_id(input$map2_shape_click$id)
  })
  observeEvent(input$map1_marker_click, {
    clicked_id(input$map1_marker_click$id)
  })
  observeEvent(input$map2_marker_click, {
    clicked_id(input$map2_marker_click$id)
  })
  
  # Update the assessment_unit dropdown when a new clicked ID is available
  observe({
    if (!is.null(clicked_id())) {
      updateSelectInput(session, 'assessment_unit', selected = clicked_id())
    }
  })
  
  # Observe changes in the assessment_unit dropdown and update clicked_id()
  observeEvent(input$assessment_unit, {
    clicked_id(input$assessment_unit)
  })
  
  # Create a reactive object for filtered_data
  df_dset_range_lf_aid <- reactive({
    if (!is.null(clicked_id()) &&
        clicked_id() != "Please make a selection" &&
        !is.null(df_dset_range_lf())) {
        debug_msg("Create a reactive object for the specific assessment unit")
      
      df_dset_range_lf() %>% 
        filter(assess_id == clicked_id())
      
    } else {
      NULL
    }
  })
  
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
  generate_ts_plot <- function(x, y, lf_select){
    
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
                       text_string = plot_label)
    
    # Convert ggplot to plotly with ggplotly
    ts_plot_plotly <- ggplotly(ts_plot, highlight = "unique_id", selected = list(marker = list(size = 10)))
    
    return(ts_plot_plotly)
  }
  
  # Generate ts_plot1
  output$ts_plot1 <- renderPlotly({
    
    if(!is.null(df_dset_range_lf_aid()) &
       !is.null(df_ken()) &
       !is.null(lifeforms())){
      debug_msg("Generate ts_plot1")
      
      ts_plot_plotly <- generate_ts_plot(x = df_dset_range_lf_aid(),
                                         y = df_ken(),
                                         lf_select = lifeforms()[1]
      )
      
      return(ts_plot_plotly)
    } else {
      return(NULL)
    }
  })
  
  # Generate ts_plot2
  output$ts_plot2 <- renderPlotly({
    
    if(!is.null(df_dset_range_lf_aid()) &
       !is.null(df_ken()) &
       !is.null(lifeforms())){
      debug_msg("Generate ts_plot2")
      
      ts_plot_plotly <- generate_ts_plot(x = df_dset_range_lf_aid(),
                                         y = df_ken(),
                                         lf_select = lifeforms()[2]
      )
      
      return(ts_plot_plotly)
    } else {
      return(NULL)
    }
  })

  # Generate a histogram of sample frequency for the time-series
  output$histo_plot <- renderPlotly({
    
    if(!is.null(df_dset_range_lf_aid())){
      debug_msg("Generate histo_plot")
      
      temp <- df_dset_range_lf_aid() 
      temp <- temp %>%
        filter(lifeform == unique(temp$lifeform)[1])
      
      histo_plot <- ggplot(temp, aes(date, n)) +
        geom_col(fill="blue") +
        scale_x_date(limits = c(min(temp$date), max(temp$date))) +
        scale_y_continuous(name="samples/month") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank()) 
      
      histo_plot_plotly <- ggplotly(histo_plot)
      return(histo_plot_plotly)
      
    } else {
      return(NULL)
    }
    
  })
  
  # Update Assessment and Comparison sliders based on range slider and filtered data
  observeEvent(df_dset_range_lf_aid(), {
    range_vals <- input$range_slider
    
    # Get the filtered data
    filtered_data <- df_dset_range_lf_aid()
    
    if (!is.null(filtered_data) & 
        !is.null(range_vals)) {
      debug_msg("Assessment and comparison slider bounds updated")
      
      # Calculate the year range of the filtered data
      filtered_year_range <- range(filtered_data$year)  # Replace 'year' with the actual variable name
      
      # Check if the filtered year range is invalid (contains -Inf and Inf)
      if (any(!is.finite(filtered_year_range))) {
        filtered_year_range <- range_vals
      }
      
      # Set Assessment slider bounds
      updateSliderInput(session, "assessment_slider", 
                        min = max(range_vals[1], filtered_year_range[1]), 
                        max = min(range_vals[2], filtered_year_range[2]))
      
      # Set Comparison slider bounds
      updateSliderInput(session, "comparison_slider", 
                        min = max(range_vals[1], filtered_year_range[1]), 
                        max = min(range_vals[2], filtered_year_range[2]))
    }
  })
  
  #generate data selection for assessment period
  df_assessment <- reactive({
    if(!is.null(df_dset_range_lf_aid()) &
       !is.null(input$assessment_slider) &
       !is.null(lifeforms())){
      debug_msg("Filter assessment data based on assessment slider input")
      
      df_temp_output <- dataSelect(df_dset_range_lf_aid(), 
                                   lims=c(min(input$assessment_slider), max(input$assessment_slider)),
                                   lf=lifeforms())
      return(df_temp_output)
    } else {
      return(NULL)
    }
  })
  
  #generate data selection for assessment period
  df_comparison <- reactive({
    if(!is.null(df_dset_range_lf_aid()) &
       !is.null(input$comparison_slider) &
       !is.null(lifeforms())){
      debug_msg("Filter comparison data based on comparison slider input")
      
      df_temp_output <- dataSelect(df_dset_range_lf_aid(), 
                                   lims=c(min(input$comparison_slider), max(input$comparison_slider)),
                                   lf=lifeforms())
      return(df_temp_output)
    } else {
      return(NULL)
    }
  })
  
  #generate data selection for assessment period
  df_assessment_qc <- reactive({
    if(!is.null(df_assessment()) &
       !is.null(lifeforms())){
      debug_msg("Quality control of assessment data")
      
      df_temp_output <- qc_assessment(df_assessment(), lf=lifeforms())
      return(df_temp_output)
    } else {
      return(NULL)
    }
  })
  
  #generate the envelope for the assessment conditions
  envelope <- reactive({
    if(!is.null(df_assessment_qc()) &
       !is.null(lifeforms()) &
       !is.null(df_comparison())){
      debug_msg("Generate the assessment envelope")
      
      #generate the qc steps for assessment data
      temp <- df_assessment_qc()
      
      #quality control steps for assessment output message
      temp_qc <- qc_output(temp)
      
      #generate the qc steps for the comparison data
      temp_comp <- qc_comparison(df_comparison())
      
      #quality control steps for the comparison output message
      temp_comp_qc <- qc_output(temp_comp)
      
      #combine the assessment and comparison qc info
      temp_combined_qc <- rbind(temp_qc, temp_comp_qc)
      
      temp <- temp %>%
        dplyr::select(assess_id, all_of(lifeforms()))
      
      #command to skip envelope fitting for data with no variance
      abort <- ifelse(
        temp_qc$tf[temp_qc$condition == "n_months_min_sufficient"]==FALSE |
          sd(as.vector(unlist(temp[,2]))) == 0 |
          sd(as.vector(unlist(temp[,3])))==0 |
          all(is.na(as.vector(unlist(temp[,2])))) |
          all(is.na(as.vector(unlist(temp[,3])))), TRUE, FALSE)
      
      if(abort==FALSE){
        
        envPts <- findEvn(as.vector(unlist(temp[,2])),
                          as.vector(unlist(temp[,3])),
                          p=0.9,
                          sc=TRUE)
        envPts_unlist <- rbindlist(envPts, fill=TRUE)
        temp_outer <- data.frame(outX=envPts_unlist$outX[complete.cases(envPts_unlist$outX)],
                                 outY=envPts_unlist$outY[complete.cases(envPts_unlist$outY)]
        )
        temp_inner <- data.frame(inX=envPts_unlist$inX[complete.cases(envPts_unlist$inX)],
                                 inY=envPts_unlist$inY[complete.cases(envPts_unlist$inY)]
        )
      } else {
        temp_outer <- NULL
        temp_inner <- NULL
      }
      
      output <- list(temp_outer, temp_inner, temp_combined_qc)
      names(output) <- c("EnvOuter", "EnvInner", "EnvQC")
      return(output)
      
    } else {
      return(NULL)
    }
  })
  
  df_pi <- reactive({
    if(!is.null(envelope()) &
       !is.null(df_comparison()) &
       !is.null(lifeforms())){
      debug_msg("Calculate the PI statistic")
      
      #use qc output from envelope function to determine whether PI can be calculated
      temp_qc <- envelope()[["EnvQC"]]
      abort <- temp_qc$tf[temp_qc$condition == "n_months_min_sufficient"]==FALSE
      
      if(abort == FALSE){
        
        temp <- df_comparison()
        compDat <- temp %>%
          dplyr::select(all_of(lifeforms())) %>%
          rename(y1=1, y2=2)
        
        env <- lapply(envelope(), function(x) x[1:2])
        piPts <- PIcalc(compDat, env, 0.9)
        piPts <- do.call(cbind.data.frame, piPts)
        colnames(piPts) <- gsub(" ", "_", colnames(piPts))
        
        print(piPts)  # Add a print statement here to check the value of df_pi()
      } else {
        
        piPts <- "The assessment data are insufficient to calculate PI statistics."
      }
      
      return(piPts)
    } else {
      return(NULL)
    }
  })
  
  # Render the text output if there is insufficient data to run the PI
  output$env_text <- renderText({
    debug_msg("Generate PI diagnostic text")
    
    if (!is.null(df_pi()) &
        class(df_pi()) == "data.frame") {
      return(NULL)  # Return NULL when df_pi() is a data frame to avoid rendering text
    } else {
      return(df_pi())  # Render the character string when df_pi() is not a data frame
    }
  })
  
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
  
  # generate the PI plot
  output$env_plot <- renderPlotly({
    if (!is.null(df_pi()) &
        !is.null(df_comparison()) &
        !is.null(clicked_id()) &
        class(df_pi()) == "data.frame") {
      debug_msg("Generate the PI plot")
      
      df_outer <- envelope()[["EnvOuter"]]
      df_inner <- envelope()[["EnvInner"]]
      names(df_outer)[1:2] <- c("x", "y")
      names(df_inner)[1:2] <- c("x", "y")
      df_outer$subid <- 1
      df_inner$subid <- 2
      df_polys <- rbind(df_inner, df_outer)
      df_polys$assess_id <- clicked_id()
      
      df_comp <- df_comparison()
      df_comp <- df_comp %>%
        rename(vx = lifeforms()[1],
               vy = lifeforms()[2]) %>%
        mutate(vx = round(vx, 2),
               vy = round(vy, 2))
      
      # grouping factor for colouring months
      df_comp$month <- as.numeric(df_comp$month)
      df_comp$season <- ifelse(df_comp$month %in% c(1, 2, 12), "12, 1, 2",
                               ifelse(df_comp$month %in% c(3, 4, 5), "3, 4, 5",
                                      ifelse(df_comp$month %in% c(6, 7, 8), "6, 7, 8",
                                             ifelse(df_comp$month %in% c(9, 10, 11), "9, 10, 11", "ERROR"))))
      
      # Create a custom color scale
      factor_levels <- unique(df_comp$season)
      myColors <- c("blue", "green", "yellow", "red")
      names(myColors) <- levels(factor_levels)
      
      pi_plot <- ggplot() +
        geom_polygon(data = subset(df_polys, subid==1), aes(x, y), fill = "grey90", colour = "black", linewidth = 0.1) +
        geom_polygon(data = subset(df_polys, subid==2), aes(x, y), fill = "white", colour = "black", linewidth = 0.1) +
        geom_path(data = df_comp, aes(x = vx, y = vy), colour = "grey", linetype = 2, linewidth = 0.25) +
        geom_point(data = df_comp, aes(x = vx, y = vy, fill = season), shape = 21, stroke=0.2, size=2) +
        scale_fill_manual(values = myColors, name = "Month") +  # Specify the fill colors and set the name for the legend
        scale_x_continuous(expand = c(0.1, 0), name = paste0("log10(",lifeforms()[1], " ", df_comp$abundance_type_units[1], ")")) +
        scale_y_continuous(expand = c(0.1, 0), name = paste0("log10(",lifeforms()[2], " ", df_comp$abundance_type_units[1], ")")) +
        ggtitle(generate_pi_label(df_pi(), df_assessment_qc(), assess_id_label=clicked_id())) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5,
                                        size = 10)) 
      
      pi_plotly <- ggplotly(pi_plot, highlight = "unique_id", selected = list(marker = list(size = 10)))
      
      return(pi_plotly)
      
    } else {
      
      return(NULL)
    }
  })
  
  # user message text to indicate the reliability of the PI
  user_prompt_text <- reactive({
    if(is.null(clicked_id())){
      output <- "Please select a map unit to generate plots or make a selection from the dropdown menu"
    } else if (!is.null(clicked_id()) &
               !is.null(envelope())) {
      
      # generate warning text for assessment period with insufficient criteria for robust assessment
      temp_qc <- envelope()[["EnvQC"]]
      temp_qc <- temp_qc %>%
        mutate(condition_long = recode(condition,
                                       n_months_min_sufficient = "The number of months in the assessment period does not fit the minimal criteria to generate the assessment envelope (i.e. n < ",
                                       n_months_robust_sufficient = "The number of months in the assessment period does not fit the robust criteria (i.e. n < ", 
                                       n_years_sufficient = "The number of years in the assessment period does not fit the robust criteria (i.e. n < ", 
                                       n_months_rep_sufficient = "The number of times each month is represented in the assessment period does not fit the robust criteria (i.e. n < ",
                                       lf_count_min_nonzero = "The proportion of zero values in the assessment period exceeds the allowable threshold for at least one lifeform (i.e. prop nonzero > ",
                                       prop_comp_monitored_min = "The proportion of unmonitored months in the comparison period exceeds the criteria for a robust assessment (i.e. prop unmonitored > "
        )
        ) %>%
        mutate(condition_long = paste0(condition_long, criterion, ").")) %>%
        filter(tf==FALSE)
      
      # if else statement to determine the initial sentence of the message displayed to the user
      if(nrow(temp_qc) > 0){
        
        if("n_months_min_sufficient" %in% temp_qc$condition){
          
          # generate first part of the message to the user
          part1 <- "The assessment data are insufficient to calculate PI statistics."
          
        } else {
          
          # generate first part of the message to the user
          part1 <- "The assessment data are sufficient to calculate PI statistics, but with low confidence due to a poorly defined assessment envelope."
          
        }
        
        # generate second part of the message to the user
        part2 <- temp_qc %>%
          summarise(conditions_comb = paste(condition_long, collapse='<br/>'))
        part2 <- unlist(part2)
        
      } else {
        
        # generate first part of the message to the user
        part1 <- "The assessment data are sufficient to present high confidence PI statistics."
        
        # generate second part of the message to the user
        part2 <- ""
        
      }
      
      # combine first and second part of the message to the user
      output_text <- paste(part1, part2, sep='<br/>')
      formatted_text <- paste0('<span style="color: black; font-size: 14px;">', output_text, '</span>')
      output <- HTML(formatted_text)
      
    } else {
      
      formatted_text <- paste0('<span style="color: red; font-size: 16px;">', 'Waiting for user input', '</span>')
      output <- HTML(formatted_text)
      
    }
    
    return(output)
  })
  
  output$user_prompt <- renderUI({
    user_prompt_text()
  })
  
  
  output$test <- renderTable({
    
    output <- head(df_pi())
    
    return(output)
  })
  
  
}
      
      
#run the app
shinyApp(ui = ui, server = server)
