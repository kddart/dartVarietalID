
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dartVarietalID

<!-- badges: start -->
<!-- badges: end -->

The goal of dartVarietalID is to …

## Installation

You can install the development version of dartVarietalID like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

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
#> Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
#> loading 'dartVarietalID'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'dartVarietalID'
#> Warning: replacing previous import 'semantic.dashboard::column' by
#> 'shiny::column' when loading 'dartVarietalID'
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
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
