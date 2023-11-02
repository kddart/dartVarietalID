#' dartVarietalIDShiny
#' @import shiny
#' @import shinyjs
#' @import shiny.semantic
#' @import semantic.dashboard
#' @import plotly
#' @import shinyWidgets
#' @import tableHTML
#' @import DT
#' @examples
#' browseURL(system.file("extdata",package = "dartVarietalID")) # Open folder of example files
#' dartVarietalIDShiny()
#' @export

dartVarietalIDShiny <- function(...) {

  structure_colors <-  c(
    "#F6222E",
    "#16FF32",
    "#FEAF16",
    "#1CFFCE",
    "#F8A19F",
    "#1C8356",
    "#85660D",
    "#BDCDFF",
    "#822E1C",
    "#3B00FB"
  )

  #### dataset list ###
  datasetListUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns('dataList'))
  }

  datasetListServer <- function(id,
                                md_list) {
    moduleServer(id,
                 function(input, output, session) {
                   output$dataList = renderUI({
                     selectizeInput(
                       inputId =  'dataList',
                       label = 'Select a sample',
                       choices =  md_list,
                       options = list(placeholder = '')
                     )
                   })
                 })
  }

  # maximum size of file to upload
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)

  box_color <- "green"
  box_height <- '100%'
  input_class <- c("ui small icon input",
                   "ui fluid icon input")

  #######################################################################
  ############### UI ####################################################
  #######################################################################

  ui <-
    dashboardPage(
      title = "Varietal Identification",

      dashboardHeader(
        color = "green",
        menu_button_label = "",
        class = "ui top attached header",
        button(
          input_id = "close",
          label = span(icon("close"), "Exit"),
          class = c("tiny",
                    "ui orange button",
                    "compact ui button")
        ),
      ),

      ### Sidebar content ###
      dashboardSidebar(
        size = "thin",
        color = "green",
        sidebarMenu(
          menuItem(
            text = span(icon("upload"), "Inputs"),
                   tabName = "inputs_tab"),
          menuItem(
            text = span(icon("map"), "References Check"),
            tabName = "ref_check"
          ),
          menuItem(
            text = span(icon("filter"), "Reference Identification"),
            tabName = "ref_id"
          ),
          menuItem(
            text = span(icon("eye"), "Visualisation"),
            tabName = "visualization"
          )
        )
      ),

      ## Body content
      dashboardBody(
        useShinyjs(),
        extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }",
                      functions = c("closeWindow")),

        #######################################################################
        ############### INPUTS ################################################
        #######################################################################

        tabItems(
          tabItem(
            tabName = "inputs_tab",
            box(
              title =  "Load/select datasets",
              color = box_color,
              width = 16,
              collapsible = FALSE,
              title_side = "top left",
              style = box_height,
              shiny.semantic::fileInput(
                inputId = "counts_file",
                label = "Load DArT Counts report",
                buttonLabel = "Browse",
                type = input_class,
                placeholder = "SEQ_SNPs_counts_0_Target.csv",
                accept = "csv"
              ),
              shiny.semantic::fileInput(
                inputId = "info_file",
                label = "Load Info file",
                buttonLabel = "Browse",
                type = input_class,
                placeholder = "InfoFile.csv",
                accept = "csv"
              ),
              shiny.semantic::fileInput(
                inputId = "snp_file",
                label = "Load DArT SNP report",
                buttonLabel = "Browse",
                type = input_class,
                placeholder = "SEQ_SNPs_0_Target.csv",
                accept = "csv"
              )
            )
          ),

          #######################################################################
          ############### REFERENCE CHECK #######################################
          #######################################################################

          tabItem(
            tabName = "ref_check",

            box(
              title = "Reference Check",
              width = 16,
              color = box_color,
              collapsible = FALSE,
              title_side = "top left",
              button(
                input_id = "run_check",
                label = span(icon("play"), "RUN"),
                class = "ui green button"
              ),
              fluidRow(
                id = "referenceValidationBody",
                title = NULL,
                hidden = TRUE,
                fluidRow(h1("Hamming Dendrogram of full dataset")),
                fluidRow(id = "referenceValidationHammingPlot",
                         column(
                           width = 12,
                           panel(style = "overflow-y:scroll; position:relative; align: centre",
                                 uiOutput("hammingPlotUi"))
                         ),
                         column(
                           width = 3,
                           textInput("searchHammingString", label = "Filter on: ", value = "")
                         )),
                fluidRow(h1("Problematic References for Review")),
                fluidRow(id = "referenceValidationParent", fluidRow(id = "referenceValidation"))
              )
            )
          ),

          #######################################################################
          ############### REFERENCE IDENTIFICATION ##############################
          #######################################################################

          tabItem(
            tabName = "ref_id",

            box(
              title = "Varietal Identification",
              width = 16,
              color = box_color,
              collapsible = FALSE,
              title_side = "top left",
              fluidRow(
                column(3,
              numericInput(
                inputId = "pop_size",
                label = "Population sizes to simulate",
                value = 10,
                min = 1
              ),
              br(),
              button(
                input_id = "run_id",
                label = span(icon("play"), "RUN"),
                class = "ui green button"
              ),
              br(),
              downloadButton('download',"Save results table to csv file"),
              br())),
              DT::dataTableOutput("res.ID"),
              style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
            )
          ),

        #######################################################################
        ############### VISUALIZATION #########################################
        #######################################################################

          tabItem(
            tabName = "visualization",

            box(
              title = "Principal Component Analysis (PCA)",
              width = 16,
              color = box_color,
              collapsible = FALSE,
              title_side = "top left",
              datasetListUI(id = "dataList"),
              button(
                input_id = "run_pca",
                label = span(icon("play"), "RUN"),
                class = "ui green button"
              ),
              fluidRow(h3("PCA of the ten closest references to the sample")),
              plotlyOutput("plot_pca")
            )
          )
        )
      )
    )

  #######################################################################
  ############### SERVER ################################################
  #######################################################################

  server <- function(input, output, session) {

    hammingCoReactive <- reactiveVal(NULL)
    ID_res <- reactiveVal(NULL)

    ### Close button ###
    observeEvent(input$close, {
      lapply(names(resourcePaths()), removeResourcePath)
      js$closeWindow()
      stopApp()
    })

    #######################################################################
    ############### REFERENCE CHECK #######################################
    #######################################################################

    observeEvent(input$run_check, {
      removeUI("#referenceValidation",
               multiple = TRUE,
               immediate = TRUE)
      insertUI(
        "#referenceValidationParent",
        where = "afterEnd",
        ui =  fluidRow(id = "referenceValidation")
      )

      infoFile = readTargetInfoFile(input$info_file$datapath)
      references = infoFile$getReferences()

      snpFilepath = input$snp_file$datapath
      snpReport = ds14.read(snpFilepath)
      snpGenotypes = ds14.genotypic(snpReport)[, rownames(references)]
      distances = dartVarietalID::dartDistanceMatrix(genotypic_data = snpGenotypes,
                                                     sampleWiseAnalysis = TRUE)

      refType_mapping = references$RefType
      names(refType_mapping) = rownames(references)
      unique_refTypes = unique(refType_mapping)

      non_distinct_refTypes_arr = lapply(unique_refTypes, function(refType) {
        targets_of_refType = names(refType_mapping)[refType_mapping %in% refType]

        max_distance_inner_group = max(distances[targets_of_refType, targets_of_refType])
        distances_between_other_refs = distances[!rownames(distances) %in%
                                                   targets_of_refType, targets_of_refType]
        min_distance_other_refs = min(distances_between_other_refs)

        if (max_distance_inner_group >= min_distance_other_refs) {
          distances_between_other_refs_invalid_only = distances_between_other_refs
          distances_between_other_refs_invalid_only[distances_between_other_refs_invalid_only >
                                                      max_distance_inner_group] = NA

          filt = apply(distances_between_other_refs_invalid_only, 1, function(x)
            any(!is.na(x)))
          references_in_violation = unique(refType_mapping[rownames(distances_between_other_refs_invalid_only)][filt])

          distances_between_other_refs_invalid_only = distances_between_other_refs_invalid_only[refType_mapping[rownames(distances_between_other_refs_invalid_only)] %in%
                                                                                                  references_in_violation, , drop = FALSE]

          non_distinct_matrix_check = distances[, targets_of_refType] <= max_distance_inner_group
          return(
            list(
              refType = refType,
              max_distance_inner_group = max_distance_inner_group,
              distances_between_other_refs_invalid_only = distances_between_other_refs_invalid_only,
              merges = sort(unique(refType_mapping[rownames(non_distinct_matrix_check)[apply(non_distinct_matrix_check, 1, function(x)
                any(x))]]))
            )
          )
        } else{
          NA
        }
      })
      non_distinct_refTypes = non_distinct_refTypes_arr[!is.na(non_distinct_refTypes_arr)]

      if (length(non_distinct_refTypes) == 0) {
        showModal(
          modalDialog(
            title = "References check",
            "All references are determined to be distinct via the hamming distance method"
          )
        )
      } else{
        insertUI("#referenceValidation",
                 where = "beforeBegin",
                 ui =
                   fluidRow(column(
                     width = 12,
                     p(
                       "Some references should be reviewed before continuing...",
                       style = "color:red"
                     ),
                     p(
                       "Tables display distances which violate the maximum genetic distance within the variety (NA indicates no violation). The row names indicate targets of the referenced variety and columns represent targets of other varities"
                     )
                   )))

        cluster_method = "complete"
        coloring =  as.character(infoFile$getReferences()$RefType)
        names(coloring) = rownames(infoFile$getReferences())
        mapping = paste(as.character(rownames(infoFile$getReferences())),
                        as.character(infoFile$getReferences()$RefType))
        names(mapping) = as.character(rownames(infoFile$getReferences()))
        hammingCo = getMergedBinsClusters(distances,
                                          cluster_method,
                                          coloring = coloring,
                                          mapping = mapping)

        hammingCoReactive(hammingCo)

        insertUI("#referenceValidation",
                 where = "afterEnd",
                 ui = div(lapply(non_distinct_refTypes, function(x) {
                   df = x$distances_between_other_refs_invalid_only
                   rownames(df) = refType_mapping[rownames(df)]
                   colnames(df) = refType_mapping[colnames(df)]

                   min_distance_for_other_var = sapply(split(df, rownames(df)), function(xx)
                     min(xx, na.rm = TRUE))

                   filterOn = unique(c(
                     rownames(x$distances_between_other_refs_invalid_only),
                     colnames(x$distances_between_other_refs_invalid_only)
                   ))

                   coloring = rep(2, length(refType_mapping))
                   coloring[refType_mapping %in% x$refType] = 1
                   names(coloring) = names(refType_mapping)
                   mapping = paste(refType_mapping, names(refType_mapping), sep = " ")
                   names(mapping) = names(refType_mapping)
                   cut = distances[rownames(distances) %in% filterOn,
                                   colnames(distances) %in% filterOn]
                   dend = getMergedBinsClusters(
                     cut,
                     cluster_method,
                     coloring = coloring,
                     mapping = mapping,
                     colorFunc = function(n)
                       c("#ff0f0f", "#878484")
                   )

                   return(fluidRow(column(
                     width = 12,
                     h2(x$refType),
                     p(
                       paste(
                         "Maximum genetic distance within the variety is ",
                         x$max_distance_inner_group
                       )
                     ),
                     renderPlot(expr = {
                       par(mar = c(12, 5, 1, 1))
                       plot(dend$dend)
                     }),
                     panel(style = "overflow-y:scroll; position:relative; align: centre", tableHTML(
                       t(x$distances_between_other_refs_invalid_only)
                     ), )
                   )))
                 })))
        show("referenceValidationBody")
      }
    })

    output$hammingPlotUi <- renderUI({
      plotOutput("hammingPlot", height = 400, width = getHammingPlotWidth())
    })

    getHammingPlotWidth <- reactive(500 + 5 * targetCount())

    targetCount <- reactive({
      hammingCo = hammingCoReactive()
      if (!is.null(hammingCo)) {
        return(2 * length(labels(hammingCo$dend)))
      } else{
        return(0)
      }
    })

    output$hammingPlot <- renderPlot({
      hammingCo = hammingCoReactive()
      if (!is.null(hammingCo)) {
        searchStr = tolower(input$searchHammingString)
        if (!is.null(searchStr) && nchar(searchStr) > 1) {
          labels_colors(hammingCo$dend) = rep("#696969", length(labels_colors(hammingCo$dend)))
          labels_colors(hammingCo$dend)[grep(searchStr, tolower(names(
            labels_colors(hammingCo$dend)
          )))] = "#ff0000"
        }
        par(mar = c(12, 5, 1, 1))
        plot(hammingCo$dend)
      } else{
        NULL
      }
    })

    #######################################################################
    ############### REFERENCE IDENTIFICATION ##############################
    #######################################################################

    observeEvent(input$run_id, {
      counts_path <- input$counts_file
      info_path <- input$info_file
      res_ID <<- runSampleAnalysis(counts.file = counts_path$datapath,
                          info.file = info_path$datapath,
                          pop.size = input$pop_size)
      ID_res(res_ID)
      output$res.ID <-  DT::renderDataTable({datatable(res_ID$res_summary)})

    })

    output$download <- downloadHandler(
      filename = function(){"dartID_results.csv"},
      content = function(fname){
        write.csv(res_ID$res_summary, fname)
      }
    )

    #######################################################################
    ############### VISUALIZATION #########################################
    #######################################################################

    observeEvent(ID_res(),{
      ID.res <- ID_res()
      output$dataList <-
        datasetListServer(id = "dataList",
                          md_list = ID.res$res_summary$TargetID.sample)
    })

    observeEvent(input$run_pca, {
      mydata_tmp <- as.character(input$dataList)
      mydata_tmp2 <-
        res_ID$res_full[which(names(res_ID$res_full) == mydata_tmp)]
      mydata_tmp3 <- mydata_tmp2[[1]][1:10, "RefType"]
      pop_ref_tmp <- res_ID$gl.references
      pop_ref <- gl.keep.pop(pop_ref_tmp,
                             pop.list = mydata_tmp3)
      pop_sam_tmp <- res_ID$gl.samples
      pop_sam <- gl.keep.pop(pop_sam_tmp,
                             pop.list = mydata_tmp)
      pop_ref_sam <- rbind(pop_ref, pop_sam)
      pop_ref_sam$other$loc.metrics <- pop_ref$other$loc.metrics

      output$plot_pca <- renderPlotly({
        pca <- gl.pcoa(pop_ref_sam, plot.out = FALSE)

        res <- gl.pcoa.plot(
          glPca = pca,
          x = pop_ref_sam,
          zaxis = 3,
          pt.size = 4,
          pt.colors = c(structure_colors[1:nPop(pop_ref)], "black"),
          axis.label.size = 2,
          label.size = 1
        )
        return(res)
      })
    })
  }

  shinyApp(ui, server, ...)
}
