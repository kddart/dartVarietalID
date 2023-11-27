
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <a href="https://www.diversityarrays.com/"><img src='www/DArT_logo.png' align="right" height="100" /></a>

# dartVarietalID

## Installation

dartVarietalID is at the moment in a private repository. To install it,
use auth_token with a token from <https://github.com/settings/tokens>.

library(devtools)

install_github(“kddart/dartVarietalID”, auth_token = “abc”)

For dartVarietalID to work the R package dartR needs to be installed.
Please consult [this installation
tutorial](https://github.com/green-striped-gecko/dartR/wiki/Installation-tutorial)
to install dartR.

## Example

This is how you can run the dartVarietalID app:

``` r
library(dartVarietalID)
dartVarietalIDShiny()
```

You can open a folder with example datasets here:

``` r
browseURL(system.file("extdata",package = "dartVarietalID"))
```

## Input Tab

<p align="center">
<img src='www/load_page.png' width="800"/>
</p>

## References Check Tab

<p align="center">
<img src='www/references_check_1.png' width="800"/>
</p>
<p align="center">
<img src='www/lreferences_check_2.png' width="800"/>
</p>

## Reference Identification Tab

<p align="center">
<img src='www/References_ID.png' width="800"/>
</p>

## Visualisation Tab

<p align="center">
<img src='www/Visualisation_1.png' width="800"/>
</p>
<p align="center">
<img src='www/Visualisation_2.png' width="800"/>
</p>
