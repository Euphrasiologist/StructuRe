# StructuRe

### Overview

A tiny package to parse STRUCTURE data to obtain the Q matrix, and then to plot the results.

### Installation and usage

Simply install using devtools (`install.packages("devtools");library(devtools)`):\
`devtools::install_github("Euphrasiologist/StructuRe")`

First, parse the data: `parsed_struct <- parseStructure("/path/to/file/")`\
Then, plot the data: `plot(parsed_struct)`\

### Acknowledgements

I acknowledge https://github.com/nicholasmason/parseStructure for the inspiration and have modified their code.
