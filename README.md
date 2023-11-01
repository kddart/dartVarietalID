
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dartVarietalID

The goal of dartVarietalID is to …

## Installation

dartVarietalID is at the moment in a private repository. To install it,
use auth_token with a token from <https://github.com/settings/tokens>.

library(devtools) install_github(“kddart/dartVarietalID”, auth_token =
“abc”)

For dartVarietalID to work the R package dartR needs to be installed.
Please consult [this installation
tutorial](https://github.com/green-striped-gecko/dartR/wiki/Installation-tutorial)
to install dartR.

## Example

This is how you can run the dartVarietalID app:

library(dartVarietalID) dartVarietalIDShiny()

You can open a folder with example datasets here:

browseURL(system.file(“extdata”,package = “dartVarietalID”))
