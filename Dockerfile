FROM r-base:latest

RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev && R -e "install.packages(c('tidyverse', 'dplyr', 'reshape2', 'growthcurver', 'lattice', 'deSolve', 'growthrates', 'ggplot2', 'ggtext', 'Hmisc'), repos='http://cran.rstudio.com/')"

RUN mkdir /workdir
WORKDIR /workdir

COPY liquid_growth_analysis.R /workdir/

ENTRYPOINT ["Rscript", "liquid_growth_analysis.R"]