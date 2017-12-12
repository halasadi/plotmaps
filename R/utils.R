
#' @title use the default eems colors
#' @useDynLib plotmaps
#' 
#' @return eems.colors vector of strings where each string represents a color
#' @export
## We need to include the useDynLib somewhere in our documentation
## see: https://stackoverflow.com/questions/36605955/c-function-not-available
default_eems_colors <- function( ) {
  writeLines(paste0("Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.\n",
                    "It combines two color schemes from the 'dichromat' package, which itself is based on\n",
                    "a collection of color schemes for scientific data graphics:\n",
                    "\tLight A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data\n",
                    "\tGraphics. EOS Transactions of the American Geophysical Union, 85(40), 385.\n",
                    "See also http://geog.uoregon.edu/datagraphics/color_scales.htm\n\n\n"))
  ## To reproduce the default eems colors:
  ## Oranges <- dichromat::colorschemes$BluetoDarkOrange.12[12:7]
  ## Blues <- dichromat::colorschemes$BrowntoBlue.12[7:12]
  ## eems.colors <- c(Oranges, "#FBFBFB", Blues)
  eems.colors <- c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", ## orange sequence
                   "#FBFBFB", ## very slightly off-white
                   "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99") ## blue sequence
  return (eems.colors)
}

