FROM rocker/shiny:4.4

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zip \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('shiny','bslib','DT','ggplot2','plotly','dplyr','tidyr','readr','jsonlite','yaml','shinyWidgets'), repos='https://cloud.r-project.org')"

WORKDIR /srv/shiny-server/HEKBlueR
COPY . /srv/shiny-server/HEKBlueR

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]

