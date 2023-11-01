
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dartVarietalID

The goal of dartVarietalID is to …

## Installation

dartVarietalID is at the moment in a private repository. To install it,
use auth_token with a token from <https://github.com/settings/tokens>.

``` r
library(devtools)
install_github("kddart/dartVarietalID", auth_token = "abc")
```

## Example

This is how you can run the dartVarietalID app:

``` r
library(dartVarietalID)
#> Loading required package: adegenet
#> Loading required package: ade4
#> 
#>    /// adegenet 2.1.10 is loaded ////////////
#> 
#>    > overview: '?adegenet'
#>    > tutorials/doc/questions: 'adegenetWeb()' 
#>    > bug reports/feature requests: adegenetIssues()
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: dartR
#> Loading required package: ggplot2
#> Loading required package: dartR.data
#> Registered S3 method overwritten by 'pegas':
#>   method      from
#>   print.amova ade4
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
#> Registered S3 method overwritten by 'genetics':
#>   method      from 
#>   [.haplotype pegas
#> Please note that rgdal will be retired during October 2023,
#> plan transition to sf/stars/terra functions using GDAL and PROJ
#> at your earliest convenience.
#> See https://r-spatial.org/r/2023/05/15/evolution4.html and https://github.com/r-spatial/evolution
#> rgdal: version: 1.6-7, (SVN revision 1203)
#> Geospatial Data Abstraction Library extensions to R successfully loaded
#> Loaded GDAL runtime: GDAL 3.5.3, released 2022/10/21
#> Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/rgdal/gdal
#>  GDAL does not use iconv for recoding strings.
#> GDAL binary built with GEOS: TRUE 
#> Loaded PROJ runtime: Rel. 9.1.0, September 1st, 2022, [PJ_VERSION: 910]
#> Path to PROJ shared files: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/rgdal/proj
#> PROJ CDN enabled: FALSE
#> Linking to sp version:1.6-1
#> To mute warnings of possible GDAL/OSR exportToProj4() degradation,
#> use options("rgdal_show_exportToProj4_warnings"="none") before loading sp or rgdal.
#> **** Welcome to dartR [Version 2.9.7 ] ****
#> Be aware that owing to CRAN requirements and compatibility reasons not all functions of the package may run after the basic installation, as some packages could still be missing. Hence for a most enjoyable experience we recommend to run the function
#> gl.install.vanilla.dartR()
#> This installs all missing and required packages for your version of dartR. In case something fails during installation please refer to this tutorial: https://github.com/green-striped-gecko/dartR/wiki/Installation-tutorial.
#> 
#> For information how to cite dartR, please use:
#> citation('dartR')
#> Global verbosity is set to: 2
#> 
#> **** Have fun using dartR! ****
#> Loading required package: parallel
#> Loading required package: data.table
#> 
#> Attaching package: 'data.table'
#> The following objects are masked from 'package:dplyr':
#> 
#>     between, first, last
#> Loading required package: Rcpp
#> Loading required package: shinyWidgets
#> Loading required package: shiny
#> Loading required package: shinyjs
#> 
#> Attaching package: 'shinyjs'
#> The following object is masked from 'package:shiny':
#> 
#>     runExample
#> The following object is masked from 'package:shinyWidgets':
#> 
#>     alert
#> The following object is masked from 'package:Rcpp':
#> 
#>     show
#> The following objects are masked from 'package:methods':
#> 
#>     removeClass, show
#> Loading required package: tableHTML
#> Loading required package: colorspace
#> Loading required package: plotly
#> 
#> Attaching package: 'plotly'
#> The following object is masked from 'package:ggplot2':
#> 
#>     last_plot
#> The following object is masked from 'package:stats':
#> 
#>     filter
#> The following object is masked from 'package:graphics':
#> 
#>     layout
#> Loading required package: semantic.dashboard
#> 
#> Attaching package: 'semantic.dashboard'
#> The following objects are masked from 'package:shiny':
#> 
#>     column, icon
#> The following object is masked from 'package:graphics':
#> 
#>     box
#> Loading required package: shiny.semantic
#> 
#> Attaching package: 'shiny.semantic'
#> The following objects are masked from 'package:semantic.dashboard':
#> 
#>     dropdown_menu, icon, menu_item
#> The following object is masked from 'package:shinyjs':
#> 
#>     toggle
#> The following objects are masked from 'package:shiny':
#> 
#>     actionButton, checkboxInput, dateInput, fileInput, flowLayout,
#>     icon, incProgress, modalDialog, numericInput, Progress,
#>     removeModal, removeNotification, selectInput, setProgress,
#>     showNotification, sliderInput, splitLayout, textAreaInput,
#>     textInput, updateActionButton, updateSelectInput,
#>     updateSliderInput, verticalLayout, withProgress
#> The following object is masked from 'package:graphics':
#> 
#>     grid
#> The following object is masked from 'package:utils':
#> 
#>     menu
#> Loading required package: DT
#> 
#> Attaching package: 'DT'
#> The following objects are masked from 'package:shiny':
#> 
#>     dataTableOutput, renderDataTable
#> Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'dartVarietalID'
#> Warning: replacing previous import 'DT::dataTableOutput' by
#> 'shiny::dataTableOutput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'semantic.dashboard::column' by
#> 'shiny::column' when loading 'dartVarietalID'
#> Warning: replacing previous import 'DT::renderDataTable' by
#> 'shiny::renderDataTable' when loading 'dartVarietalID'
#> Warning: replacing previous import 'semantic.dashboard::icon' by 'shiny::icon'
#> when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::textInput' by
#> 'shiny.semantic::textInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::showNotification' by
#> 'shiny.semantic::showNotification' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::incProgress' by
#> 'shiny.semantic::incProgress' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::modalDialog' by
#> 'shiny.semantic::modalDialog' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::removeModal' by
#> 'shiny.semantic::removeModal' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::sliderInput' by
#> 'shiny.semantic::sliderInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::selectInput' by
#> 'shiny.semantic::selectInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'semantic.dashboard::dropdown_menu' by
#> 'shiny.semantic::dropdown_menu' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::setProgress' by
#> 'shiny.semantic::setProgress' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::fileInput' by
#> 'shiny.semantic::fileInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::withProgress' by
#> 'shiny.semantic::withProgress' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::verticalLayout' by
#> 'shiny.semantic::verticalLayout' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::numericInput' by
#> 'shiny.semantic::numericInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::actionButton' by
#> 'shiny.semantic::actionButton' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::updateSliderInput' by
#> 'shiny.semantic::updateSliderInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::updateActionButton' by
#> 'shiny.semantic::updateActionButton' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::removeNotification' by
#> 'shiny.semantic::removeNotification' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::textAreaInput' by
#> 'shiny.semantic::textAreaInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::dateInput' by
#> 'shiny.semantic::dateInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'semantic.dashboard::menu_item' by
#> 'shiny.semantic::menu_item' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::Progress' by
#> 'shiny.semantic::Progress' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::checkboxInput' by
#> 'shiny.semantic::checkboxInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::splitLayout' by
#> 'shiny.semantic::splitLayout' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::icon' by 'shiny.semantic::icon' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::updateSelectInput' by
#> 'shiny.semantic::updateSelectInput' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::flowLayout' by
#> 'shiny.semantic::flowLayout' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny.semantic::toggle' by
#> 'shinyjs::toggle' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shinyWidgets::alert' by 'shinyjs::alert'
#> when loading 'dartVarietalID'
#> Warning: replacing previous import 'methods::removeClass' by
#> 'shinyjs::removeClass' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny::runExample' by 'shinyjs::runExample'
#> when loading 'dartVarietalID'
#> Warning: replacing previous import 'methods::show' by 'shinyjs::show' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'dartVarietalID'
#> Warning: replacing previous import 'plotly::filter' by 'stats::filter' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'Rcpp::.DollarNames' by
#> 'utils::.DollarNames' when loading 'dartVarietalID'
#> Warning: replacing previous import 'shiny.semantic::menu' by 'utils::menu' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'Rcpp::prompt' by 'utils::prompt' when
#> loading 'dartVarietalID'
dartVarietalIDShiny()
#> PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
```

<div style="width: 100% ; height: 400px ; text-align: center; box-sizing: border-box; -moz-box-sizing: border-box; -webkit-box-sizing: border-box;" class="muted well">Shiny applications not supported in static R Markdown documents</div>

You can open a folder with example datasets here:

``` r
browseURL(system.file("extdata",package = "dartVarietalID"))
```
