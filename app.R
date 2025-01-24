# Increase the maximum upload size to 160 MB 
options(shiny.maxRequestSize = 100*1024^2)

#install and/or load the required packages
library("data.table")
library("fst")
library("leaflet")
library("dplyr")
library("shiny")
library("sf")
library("tidyverse")
library("broom")
library("dichromat")
library("pracma")
library("shinyjs")
library("plotly")
library("shinyWidgets")
library("trend")
library("scales")
library("EnvStats")
library("stringr")

#enter the filename of the multipolygon shapefile the data is partitioned with
file_shp_part <- "COMP4_WFD.shp"

#enter the filename of the point shapefile the station data is partitioned with
file_pts_part <- "Stations.shp"

#enter the main directory to use to access the processed data
dir_raw <- "data/"

#generate a string for the directory of the multipolygon shapefile
dir_shp <- paste0(dir_raw, gsub(".shp","",file_shp_part), "/")

#generate a string for the directory of the point shapefile
dir_pts <- paste0(dir_raw, gsub(".shp","",file_pts_part), "/")

#load the supporting functions for calculating the indicator
source("pi_functions.R")

#load the supporting functions for everything else
source("app_functions.R")

#read the plankton data
df_plot <- read_fst(path=paste0(dir_raw, "COMP4_WFD_Stations_plankton", ".fst")) %>%
  filter(tconf >= 0.5) %>%
  group_by(data_id, assess_id, lifeform) %>%
  filter(length(unique(year)) >= 3) %>%
  ungroup()
  
#load the shapefile associated with the data
shp_part <- st_read(paste(dir_shp, file_shp_part, sep=""))

#generate points shape file to represent the stations
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
  df_assessment <- reactiveVal(NULL)
  df_comparison <- reactiveVal(NULL)
  df_assessment_qc <- reactiveVal(NULL)
  envelope <- reactiveVal(NULL)
  df_pi <- reactiveVal(NULL)
  
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
          dplyr::select(lifeform, assess_id, is_point, n, estimate1, p.value, sens_estimate, sens_p.value) %>%
          dplyr::rename(statistic = estimate1)
        
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
  
  # Generate ts_plot1
  output$ts_plot1 <- renderPlotly({
    
    if(!is.null(df_dset_range_lf_aid()) &
       !is.null(df_ken()) &
       !is.null(lifeforms())){
      debug_msg("Generate ts_plot1")
      
      ts_plot_plotly <- generate_ts_plot(x = df_dset_range_lf_aid(),
                                         y = df_ken(),
                                         lf_select = lifeforms()[1],
                                         assessment_range = input$assessment_slider,
                                         comparison_range = input$comparison_slider
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
                                         lf_select = lifeforms()[2],
                                         assessment_range = input$assessment_slider,
                                         comparison_range = input$comparison_slider
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
      
      suppressWarnings({
        histo_plot <- ggplot() +
          geom_col(data=temp, aes(date, n, text = paste0("Samples: ", n, "<br>",
                                                         "Month: ", format(date, "%Y-%m"))),
                   fill="blue") +
          scale_x_date(limits = c(min(temp$date), max(temp$date))) +
          scale_y_continuous(name="samples/month") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5),
                axis.title.x = element_blank())
      })
      
      histo_plot_plotly <- ggplotly(
        histo_plot,
        tooltip = "text"  # Use the "text" aesthetic for tooltips
      )
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
      filtered_range <- range(filtered_data$year)  # Replace 'year' with the actual variable name
      
      # Check if the filtered year range is invalid (contains -Inf and Inf)
      if (any(!is.finite(filtered_range))) {
        filtered_range <- range_vals
      }
      
      #calculate range values
      min_val <- max(range_vals[1], filtered_range[1]) 
      max_val <- min(range_vals[2], filtered_range[2])
      
      #logical statement for if range is over ten then assess range and comp range are 5 years
      if(max_val - min_val >= 9){
        assess_range_vals <- c(min_val, min_val+4)
        comp_range_vals <- c(max_val-4, max_val)  
      } else {
        int_val <- (max_val - min_val) %/% 2
        assess_range_vals <- c(min_val, min_val + int_val)
        comp_range_vals <- c(max_val - int_val, max_val)
      }
      assess_range_vals
      comp_range_vals

      # Set Assessment slider bounds
      updateSliderInput(session, "assessment_slider", 
                        value=assess_range_vals,
                        min=min_val,
                        max=max_val)
      
      # Set Comparison slider bounds
      updateSliderInput(session, "comparison_slider", 
                        value=comp_range_vals,
                        min=min_val,
                        max=max_val)
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
      temp_assess <- df_assessment_qc()
      
      #quality control steps for assessment output message
      temp_assess_qc <- qc_text_output(temp_assess)
      
      #generate the qc steps for the comparison data
      temp_comp <- qc_comparison(comp = df_comparison(), assess = temp_assess)
      
      #quality control steps for the comparison output message
      temp_comp_qc <- qc_text_output(temp_comp)
      
      #combine the assessment and comparison qc info
      temp_combined_qc <- rbind(temp_assess_qc, temp_comp_qc)
      
      temp_assess <- temp_assess %>%
        dplyr::select(assess_id, all_of(lifeforms()))
      
      #command to skip envelope fitting for data with no variance
      abort <- ifelse(
        temp_combined_qc$tf[temp_combined_qc$condition == "n_months_min_sufficient"]==FALSE |
          temp_combined_qc$tf[temp_combined_qc$condition == "overlap_not_present"]==FALSE |
          sd(as.vector(unlist(temp_assess[,2]))) == 0 |
          sd(as.vector(unlist(temp_assess[,3])))==0 |
          all(is.na(as.vector(unlist(temp_assess[,2])))) |
          all(is.na(as.vector(unlist(temp_assess[,3])))), TRUE, FALSE)
      
      if(abort==FALSE){
        
        envPts <- findEvn(as.vector(unlist(temp_assess[,2])),
                          as.vector(unlist(temp_assess[,3])),
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
      abort <- temp_qc$tf[temp_qc$condition == "n_months_min_sufficient"]==FALSE |
        temp_qc$tf[temp_qc$condition == "overlap_not_present"]==FALSE
      
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
        
        piPts <- "The data are insufficient or unsuitable for calculating meaningful PI statistics."
      }
      
      return(piPts)
    } else {
      return(NULL)
    }
  })
  
  # Render the text output if there is insufficient data to run the PI
  output$env_text <- renderText({
    if (!is.null(df_pi()) &
        class(df_pi()) == "data.frame") {
      debug_msg("Generate PI diagnostic text")
      
      return(NULL)  # Return NULL when df_pi() is a data frame to avoid rendering text
    } else {
      return(df_pi())  # Render the character string when df_pi() is not a data frame
    }
  })
  
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
        dplyr::rename(vx = lifeforms()[1],
               vy = lifeforms()[2]) %>%
        mutate(vx = round(vx, 2),
               vy = round(vy, 2))
      
      seasons <- c("Winter", "Spring", "Summer", "Autumn")
      
      # grouping factor for colouring months
      df_comp$month <- as.numeric(df_comp$month)
      df_comp$season <- ifelse(df_comp$month %in% c(1, 2, 12), seasons[1],
                               ifelse(df_comp$month %in% c(3, 4, 5), seasons[2],
                                      ifelse(df_comp$month %in% c(6, 7, 8), seasons[3],
                                             ifelse(df_comp$month %in% c(9, 10, 11), seasons[4], "ERROR"))))
      df_comp$season <- factor(df_comp$season, levels=seasons)
      
      # Create a custom color scale
      myColors <- c("blue", "green", "yellow", "red")
      names(myColors) <- seasons
      
      units <- df_comp$abundance_type_units[1]
      
      suppressWarnings({
        pi_plot <- ggplot() +
          geom_polygon(
            data = subset(df_polys, subid == 1),
            aes(x, y),  
            fill = "grey90", colour = "black", linewidth = 0.1
          ) +
          geom_polygon(
            data = subset(df_polys, subid == 2),
            aes(x, y),  
            fill = "white", colour = "black", linewidth = 0.1
          ) +
          geom_path(
            data = df_comp,
            aes(x = vx, y = vy),  # Customize tooltip for this layer
            colour = "grey", linetype = 2, linewidth = 0.25
          ) +
          geom_point(
            data = df_comp,
            aes(x = vx, y = vy, fill=season, text = paste0("log10(", lifeforms()[1], " ", units, ")", ": ", vx, "<br>",
                                                             "log10(", lifeforms()[2], " ", units, ")", ": ", vy, "<br>",
                                                             "Month: ", format(as.Date(paste(year, month, 15, sep="-")), "%Y-%m"))),
            stroke = 0.2, size = 2, shape = 21
          ) +
          scale_fill_manual(values = myColors, name = "Season") +
          scale_x_continuous(expand = c(0.1, 0), name = paste0("log10(", lifeforms()[1], " ", units, ")")) +
          scale_y_continuous(expand = c(0.1, 0), name = paste0("log10(", lifeforms()[2], " ", units, ")")) +
          ggtitle(generate_pi_label(df_pi(), df_assessment_qc(), assess_id_label = clicked_id())) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 10))
      })
      
      pi_plotly <- ggplotly(
        pi_plot,
        tooltip = "text"  # Use the "text" aesthetic for tooltips
      )
      
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
                                       overlap_not_present = "The assessment and comparison periods are overlapping and an accurate PI cannot be calculated (i.e. overlapping months > ",
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
        
        if("n_months_min_sufficient" %in% temp_qc$condition |
           "overlap_not_present" %in% temp_qc$condition){
          
          # generate first part of the message to the user
          part1 <- "The data are insufficient or unsuitable for calculating meaningful PI statistics."
          
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
  
  #####################################################################################################################
  # end of necessary server functions
  #####################################################################################################################
  
  #generate diagnostic output
  test <- reactive({
    
    #use qc output from envelope function to determine whether PI can be calculated
    #temp_qc <- envelope()[["EnvQC"]]
    #return(temp_qc)
    
    df_comp <- df_comparison()
    df_comp <- df_comp %>%
      dplyr::rename(vx = lifeforms()[1],
                    vy = lifeforms()[2]) %>%
      mutate(vx = round(vx, 2),
             vy = round(vy, 2))
    
    # grouping factor for colouring months
    df_comp$month <- as.integer(df_comp$month)
    df_comp$season <- ifelse(df_comp$month %in% c(1, 2, 12), "Winter",
                             ifelse(df_comp$month %in% c(3, 4, 5), "Spring",
                                    ifelse(df_comp$month %in% c(6, 7, 8), "Summer",
                                           ifelse(df_comp$month %in% c(9, 10, 11), "Autumn", "ERROR"))))
    
    return(df_comp)
    
  })

  
  output$test <- renderTable({
    
    output <- test()
    
    return(output)
  })
  
  
}
      
      
#run the app
shinyApp(ui = ui, server = server)
