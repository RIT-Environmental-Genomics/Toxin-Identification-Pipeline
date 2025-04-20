# DOCKER

For the use in this experiment, a [Jupyter Docker](https://jupyter-docker-stacks.readthedocs.io/en/latest/) environment was implemented onto an [Ubuntu 22.04](https://www.proxmox.com/en/) LXC Container in [Proxmox](https://www.proxmox.com/en/)

<p>
The following Hardware was used to drive this project. Note that the total hardware is not the same as what may be required for this project. A lower and upper baseline has not at this time been identified.
<br >
<br >

|Hardware|Information|
| ------ | ------ |
|Server|Dell PowerEdge R710 Server Rack|
|CPU|2x 12 Core Intell Xeon 2.13 GHZ|
|RAM|196 GB DDR3 2200 MHZ|
|GPU|N/A|
|Storage|6TB HDD|

  <br >
  <br >
</p>


 ## WARNING: IT IS RECOMMENDED YOU USE A NON-ROOT USER THAT YOU PROVIDE SUDO TO BUILD ALL ENVIRONMENTS
<p>
  <br>
  
Example
```sh
adduser conda
```
<br >
<br >

Following should appear:
<br >

```sh
Adding user `conda' ...
Adding new group `conda (1000) ...
Adding new user `conda' (1000) with group `conda' ...
Creating home directory `/home/conda' ...
Copying files from `/etc/skel' ...
New Password: []
```
<br >
Go ahead and add a memorable password, 
<br >
<br >
The following will now appear:
<br >
<br >

```sh
Enter the new value, or press ENTER for the default
        Full Name []: 
        Room Number []: 
        Work Phone []: 
        Home Phone []: 
        Other []: 
Is the information correct? [Y/n]
```

No information is required for this but can be put in if desired
<br >
<br >
Now give user conda the sudo privileges:
<br >
<br >
```sh
usermod -aG sudo conda
```
<br >
<br >
</p>

 ## Dockerize:

<p>
<br >
Jupyter was Dockerized through Dockerfile using the following process:
  <br >
   <br > 

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
<br >
<br >
</p>


## Environment
<p>
<br >
Docker Environment was then built using docker's build command:

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
sudo docker run -p <active port>:<listening port> <Container Name>
```
<br > 
<br >
   <br > 
Example of command:
<br > 

```sh
sudo docker run -d -p 8888:8888 jupyter
```

<br >
<br >
</p>

### Warning: additional information: 

<p>
  
When using -d to run docker as this example has, you MUST check your logs to identify the token for your Jupyter container
<br >
To do this, simplyrun docker container command:

```sh
sudo docker container ls -a
```
To show all running and not running containers
<br >
If this is your first time running Docker, there should only be one container (your "jupyter" container) available, 
<br >
<br >
Example:
<br >
```sh
CONTAINER ID   IMAGE     COMMAND                  CREATED      STATUS        PORTS                                         NAMES
3936aee05179   jupyter   "bash -c 'source /op…"   4 days ago   Up 24 hours   0.0.0.0:8888->8888/tcp, [::]:8888->8888/tcp   hardcore_allen
```

To start, stop, reset or use logs, you will use the container ID of this jupyter container. If a container ever fails it can be identified as such under the status section and can be restarted. Do no create a new container if one already exists, it will not contain previous information and will require you to reinstall all previously installed softwares. simply back up previous saves consistently and use the log command liberally.

Using log:

```sh
sudo docker container logs 3936aee05179
```

Logfiles will be provided, you will need to identify the following within the list (Example):

```sh
[I 2025-04-19 14:17:41.494 ServerApp] http://3936aee05179:8888/tree?token=cb712a3c5c689af38a240f31e5f28b16671705a8d16657fe
[I 2025-04-19 14:17:41.494 ServerApp]     http://127.0.0.1:8888/tree?token=cb712a3c5c689af38a240f31e5f28b16671705a8d16657fe
```

The provided token at the end will be your token to acces Jupyter. Please be sure when accessing jupyter you use the provided IP address of your device/network/container in place of 127.0.0.1. you will still provide the port to listen to (in this case 8888) and the remaining portions of the string. 


</p>
