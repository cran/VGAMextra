##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.


.onAttach <- function(libname = NULL, pkgname = "VGAMextra") {
  packageStartupMessage("\n     =====    VGAMextra 0.0-4    ===== \n\n",
                        "Additions and extensions of the package VGAM.",
                        "\n",
                        "For more on VGAMextra, visit" , "\n",
                "     https://www.stat.auckland.ac.nz/~vmir178/",  
                "\n\n", 
                "For a short description, fixes/bugs, and new" ,
                "\n",
                "features type vgamextraNEWS().")
  invisible()
}

vgamextraNEWS <- function() {
  file.show(system.file("doc",
                        "sceneVGAMextra.txt", package="VGAMextra"))
}
