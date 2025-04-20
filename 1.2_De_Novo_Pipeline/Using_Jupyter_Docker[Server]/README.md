# DOCKER

For the use in this experiment, a [Jupyter Docker](https://jupyter-docker-stacks.readthedocs.io/en/latest/) environment was implemented onto an [Ubuntu 22.04](https://www.proxmox.com/en/) LXC Container in [Proxmox](https://www.proxmox.com/en/)
<p>
The following Hardware was used to drive this project. Note that the total hardware is not the same as what may be required for this project. A lower and upper baseline has not at this time been identified.
<br >
  <br >
</p>

|Hardware|Information|
| ------ | ------ |
|Server|Dell PowerEdge R710 Server Rack|
|CPU|2x 12 Core Intell Xeon 2.13 GHZ|
|RAM|196 GB DDR3 2200 MHZ|
|GPU|N/A|
|Storage|6TB HDD|

<p>
<br >
  Miniconda was pulled using the following docker command:
  <br > 
   <br > 
</p>

```sh
docker pull continuumio/miniconda3
```

<p>
<br >
Jupyter was Dockerized through Dockerfile using the following process:
  <br >
   <br > 
</p>

```sh
FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && \
    apt-get install -y wget bzip2 git curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget -qO /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set PATH
ENV PATH="/opt/conda/bin:$PATH"

# Change default shell to bash for later RUN and CMD
SHELL ["/bin/bash", "-c"]

# Install packages into base environment
RUN conda update -n base -c defaults -y conda && \
    conda install -y jupyter

# Make sure Jupyter lab is working with conda, was having issues initializing conda base
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "conda activate base" >> /etc/bash.bashrc

# Standard Jupyter port, can change if you wish
EXPOSE 8888

# Set default command to run Jupyter Notebook
CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate base && jupyter notebook --ip=0.0.0.0 --allow-root"]

# IP=0.0.0.0 whitelists all IP addresses to allow anyone to join the jupyter notebook
```

