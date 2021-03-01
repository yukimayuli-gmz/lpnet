# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Print hello world
hello <- function() {
  print("Hello, world!")
}

#4.1
usethis::use_package("ape")
#ape::nj()#使用函数nj
#5.1
#devtools::load_all()
#devtools::document()
