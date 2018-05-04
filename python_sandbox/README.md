# Overview
The goal of this document is to enable the correct python modules (packages of scripts that are used by other programs) that we will need in the course.  The modules are:
1. numpy
2. scipy
3. Forthon
4. matplotlib
as well as to connect the correct version of python, 2.7.X, to the codes we will be running.  If you already know how to do this, by all means, do it however you want.  Otherwise, this is how I have accomplished this goal.
If you know 

# Creating the Sandbox (1 time only)  

I like to create my python environment with a python package called VirtualEnv. More about this program can be read [here](https://virtualenv.pypa.io/en/stable/).  Before you even start the installation process though, make sure your python build is 2.7.X (check with `~/path/python --version`).  If the output is something other than 2.7.X, where X is any number, you will need to build a version of python 2.7.X.  You can find directions on how to do this online, but it often depends on specifics of your computer's opertating systems.  If you need help, feel free to email us.

The installation procedure I use is the following:

1. I create a hidden directory called .virtual_envs in my home directory, i.e. `mkdir ~/.virtual_envs` for Linux or OSX.

2.  I go to this directory, i.e. `cd ~/.virtual_envs` for Linux or OSX.

3.  In this directory, I follow the installation commands found [here](https://virtualenv.pypa.io/en/stable/installation/).
Specifically
```bash
curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-15.1.0.tar.gz
tar xvfz virtualenv-15.1.0.tar.gz
~/path/python virtualenv-15.1.0/virtualenv.py warp
```

where `~/path` indicates the path to the python you wish to be associated with the environment.  The above commands will create the directory `~/.virtual_envs/warp` and a bunch of materials therein.

# Turning the sandbox on and off: 
  ## On 
  ```source ~/.virtual_envs/warp/bin/activate```  
  ## Off 
  ```deactivate``` 

# Setting up standard libraries (1 time only)
1.  Turn your sandbox on, `source ~/.virtual_envs/warp/bin/activate`.
2.  Pip is a package manager that is automatically installed by virtualenv.  I install my python modules with pip.  Issue the following command:

```
pip install numpy scipy Forthon matplotlib
```

Alternatively, you can install each module in succession
```
pip install numpy
pip install scipy
pip install Forthon
pip install matplotlib
```
this will install the packages numpy, scipy, Forthon, and matplotlib.  All of these packages except Forthon are standard packages used by computational scientists, and they have substantial, accessible documentation online both formally and informally on Q&A forums such as StackExchange.  Forthon, on the other hand, is a bit of a black box to me with almost no documentation and almost no community Q&A support.  It was written by the people who wrote Warp to handle the communication between the slower python code you will be working with and the much faster Fortran algorithms, and is therefore a much smaller-community package.
