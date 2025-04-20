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


 ## WARNING: IT IS RECOMMENDED YOU USE A NON-ROOT USER THAT YOU PROVIDE SUDO TO BUILD ALL ENVIRONMENTS

an example
```sh
adduser conda
```
Following should appear:

```sh
Adding user `conda' ...
Adding new group `conda (1002) ...
Adding new user `conda' (1002) with group `conda' ...
Creating home directory `/home/conda' ...
Copying files from `/etc/skel' ...
New Password: []
```
go ahead and add a memorable password

the following will now appear:

```sh
Enter the new value, or press ENTER for the default
        Full Name []: 
        Room Number []: 
        Work Phone []: 
        Home Phone []: 
        Other []: 
Is the information correct? [Y/n]
```
no information is required for this but can be put in if desired

now give user conda the sudo privileges:

```sh
usermod -aG sudo conda
```

 ## Dockerize:

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

## Environment

Docker Environment was then built using docker's build command:
<p>
  <br >

```sh
sudo docker build -t <Container Name> /path/to/Dockerfile
```
<br >
   <br > 
</p>

Example of this command used for this project:

```sh
sudo docker build -t jupyter .
```

<p>
  <br >

Finally Docker run was done with port [-p] 8888:8888 accessable to allow local (or when port forwarded, global) access to the jupyter container using a token authenticator built into Jupyter [No I wasn't going to set up a VPN that's a future Andrew's problem]

```sh
sudo docker run -d -p <active port>:<listening port> <Container Name>
```
<br >
   <br > 
</p>

Example of used commande:
```sh
sudo docker run -d -p 8888:8888 jupyter
```
