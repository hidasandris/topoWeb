# topoWeb

An R package to calculate topological importance (TI, WI), and interaction asymmetry.

## Installation

-   With devtools: `devtools::install_github("hidasandris/topoWeb")`
-   With pak: `pak::pak("hidasandris/topoWeb")`

## Usage

### TI/WI

Calculate the topological importance (TI) and it's weighted version (WI) in case of weighted network.

``` r
toy_data |>
calculate_TI_WI(3)
```

This function computes the TI index by considering 3 nodes and accounting for indirect effects (Jordán 2009). Note that the WI index is not calculated, as the toy_data lacks weight attributes.

You can select the nodes with the highest *n* TI values with:

``` r
toy_data |>
calculate_TI_WI(3) |>
dplyr::slice_max(TI_index, n = 5)
```

### Interaction asymmetry

Calculate the interaction asymmetry of a network.

``` r
toy_data |>
calculate_asymmetry(filter = F)
```

This function returns directed links and the weight of asymmetric effects (Jordán 2024). The most asymmetric interactions in a network can be a proxy of the most important regulation effects. In case of food webs the method can indicate the strongest top-down, bottom-up or indirect controls.

### References

Jordán, F., 2009. Keystone species and food webs. Philos. Trans. R. Soc. B Biol. Sci. 364, 1733–1741. <https://doi.org/10.1098/rstb.2008.0335>

Jordán, F., Capelli, G., Primicerio, R., Bodini, A., 2024. Strongly asymmetric interactions and control regimes in the Barents Sea: a topological food web analysis. Front. Mar. Sci. 11. <https://doi.org/10.3389/fmars.2024.1301612>
