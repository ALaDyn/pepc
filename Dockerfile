# lean base image to start from
FROM alpine:3.17

# add labels
LABEL description="An image to plot and publish results from JUBE runs"
LABEL version="1.1"

# install gnuplot, pandas, plotly, and mkdocs
RUN apk update && \
    apk add \
    gnuplot \
    sqlite \
    git \
    curl \
    python3 \
    jupyter-notebook \
    py3-pandas \
    py3-pip ;\
    pip3 install plotly ;\
    pip3 install pybadges mkdocs mkdocs-material mkdocs-jupyter

# try and have start dir?
WORKDIR /tmp
