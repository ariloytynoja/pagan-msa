Installation of PAGAN
=====================

![](http://wasabiapp.org/download/theme/icons/pagan_logo.png)

PAGAN is developed and most thoroughly tested on Ubuntu GNU/Linux system.

*   [Pre-compiled executables for Linux](#pre-compiled-executables-for-linux)
*   [Pre-compiled executables for Windows](#pre-compiled-executables-for-windows)
    *   [Hints for using PAGAN on Windows](#hints-for-using-pagan-on-windows)
*   [Pre-compiled executables for Mac OSX](#pre-compiled_executables_for_Mac_OSX)
*   [Using PAGAN2 with Docker](#using-pagan2-dwith-docker)
*   [Installation of PAGAN from source code](#installation-of-pagan-from-source-code)
    *   [Installation of Boost libraries](#installation-of-boost-libraries)
    *   [Download and installation of PAGAN using git](#download-and-installation-of-pagan-using-git)
    *   [Installation of Exonerate, MAFFT, RAxML](#installation-of-helper-tools)
*   [Advanced: installing boost from source code](#advanced:-installing-boost-from-source-code)

[Back to PAGAN home.](../README.md)

### Pre-compiled executables for Linux

The pre-compiled executables are available in the [binaries folder](../binaries/).

The PAGAN executable for Linux has been compiled on Red Hat ELS 6.5 (64bit) and confirmed to work also on Ubuntu 13.10 and CentOS 6.0. The package contains the PAGAN executable as well as the helper programs needed for the full functionality.

Using the command line, the file can be downloaded and unpacked with the commands:

```
mkdir ~/programs
cd ~/programs
wget http://wasabiapp.org/download/pagan/pagan.linux64.[latest].tgz
tar xvzf pagan.linux64.[latest].tgz
./pagan/bin/pagan
```

For the ease of use, it is recommended to add the directory pagan/bin to the system path or copy its content to a such directory (such as ~/bin or /usr/local/bin). This can be done with commands similar to these:

```
echo "export PATH=$PATH:/home/$USER/programs/pagan/bin" >> ~/.bash_profile
```

or

```
cp -R /home/$USER/programs/pagan/bin/* ~/bin/
```

Note that if the executables of the bundled version are moved or copied, the directory structure (```bin/lib```) and the location of library dependencies have to be retained.

### Pre-compiled executables for Windows

Currently pre-compiled binaries of the latest version of PAGAN are provided for Linux only. The pre-compiled executables of some older versions of PAGAN are available in the [binaries folder](../binaries/).

It is recommended to place the unpacked content of the zip file to the directory “Program Files”. Here we assume that the file path to the directory “pagan” is ```C:\Program Files\pagan```.

For the ease of use, it is recommended to add the directory pagan/bin permanently to the system path. On Windows 7, this can be done through the following steps:

*   open “Control Panel”
*   choose “System”
*   on the left panel, choose “Advanced system settings”
*   on the new window, click “Environment Variables”
*   find “Path” in the lower field, select it and click “Edit”
*   leave the existing content intact and add **;C:\\"Program Files"\\pagan** in the end of the field “Variable value”
*   click “OK” three times and close “Control Panel”

If you now open Command Prompt and write pagan, the program should run and show the options list.

The executable has been compiled and tested on Windows 7 (32bit) system running in a virtual machine.

#### Hints for using PAGAN on Windows

PAGAN is used through commands given in Command Prompt. For those not familiar with moving around files, the following workflow may be useful.

 *   open “Command Prompt” from “All Programs” -> “Accessories”
 *   in Command Prompt, type commands

```
mkdir analysis
cd analysis
explorer .
```

*   drag and drop your alignment input file(s) to the Explorer window that was opened
*   in Command Prompt, run PAGAN as ```pagan --seqfile=myfile [other options]```
*   once finished, the result files will appear in the Explorer window

### Pre-compiled executables for Mac OSX

Currently pre-compiled binaries of the latest version of PAGAN are provided for Linux only. The pre-compiled executables of some older versions of PAGAN are available in the [binaries folder](../binaries/).

The PAGAN executable for Mac has been compiled on OSX 10.7.5. The package contains the PAGAN executable as well as the helper programs needed for the full functionality.

The PAGAN package file can be downloaded and unpacked using graphical software. With the command line, the same can be done with the commands:

```
mkdir ~/programs
cd ~/programs
curl -O http://wasabiapp.org/download/pagan/pagan.osx64.[latest].tgz
tar xvzf pagan.osx64.[latest].tgz
./pagan/bin/pagan
```

For the ease of use, it is recommended to add the directory pagan/bin to the system path or copy its content to a such directory (such as ```~/bin``` or ```/usr/local/bin```). This can be done with commands similar to these:

```
echo "export PATH=$PATH:/Users/$USER/programs/pagan/bin" >> ~/.bash_profile
```

or

```
cp -R /Users/$USER/programs/pagan/bin/* ~/bin/
```

Note that if the executables are moved or copied, the directory structure (bin/mafftlib) and the location of library dependencies have to be retained.

### Using PAGAN2 with Docker

PAGAN2 for Linux is provided as Docker image that can be used on all computer platforms supporting [Docker](https://www.docker.com/).  
 

**Download PAGAN2 image:**

```
docker pull ariloytynoja/pagan2
docker tag ariloytynoja/pagan2 pagan2
```

**Run the image:**

```
docker run --rm -v $path_to_current_directory:/data pagan2 -s infile.fas -o outfile
```

**Alternative: Create a helper script (for Linux or OSX) and run that:**

```
cat > pagan2.sh << EOF 
#!/bin/bash
docker run --rm -v `pwd`:/data pagan2 "\$@"
EOF

chmod +x pagan2.sh 

./pagan2.sh -s infile.fas -o outfile
```

### Installation of PAGAN from source code

#### Installation of Boost libraries

PAGAN requires *-dev versions* of three utility libraries from [http://www.boost.org](http://www.boost.org) that may not be included in standard OS installations. These libraries have to be installed **before** compiling the PAGAN source code. On Ubuntu, they can be installed using commands:

```
apt-cache search libboost-program-options
apt-cache search libboost-regex
apt-cache search libboost-system
apt-cache search libboost-thread
```

(See which version number, ending with -dev, is provided and edit the command below.)

```
sudo apt-get install libboost-program-options1.54-dev
sudo apt-get install libboost-regex1.54-dev
sudo apt-get install libboost-system1.54-dev
sudo apt-get install libboost-thread1.54-dev
```

If your distribution does not provide Boost libraries (highly unlikely) or you are not allowed to install software on your system, you can download and install the necessary library from the Boost project repositories directly into your PAGAN source code directory. See the Advanced material below for details.

On MacOSX, Boost can be installed with package management systems such as [Homebrew](http://mxcl.github.com/homebrew/). On OSX 10.7, it can be installed with the following commands:

```
# uncheck below to install Homebrew
# sudo /usr/bin/ruby -e "$(curl -fsSL https://raw.github.com/gist/323731)"
sudo brew install boost
```

See [http://mxcl.github.com/homebrew/](http://mxcl.github.com/homebrew/) for further details.

#### Download and installation of PAGAN using git

The most recent version of the PAGAN source code is available from the git-repository and snapshots of this are downloadable as compressed tar-packages. The git software can be found at [http://git-scm.com](http://git-scm.com). On Ubuntu, it can also be installed using command:

```
sudo apt-get install git-core
```

Given that git is installed, the PAGAN source code can be downloaded and compiled using commands:

```
git clone https://github.com/ariloytynoja/pagan-msa.git
cd pagan-msa/src
make
./pagan
```

Without git, the latest (but possibly not up-to-date) tar-packaged source code can be downloaded and compiled using commands:


```
wget https://github.com/ariloytynoja/pagan-msa/blob/master/binaries/pagan.src.[latest].tgz
tar xzf pagan.src.[latest].tgz
cd pagan-msa/src
make
./pagan
```


#### Installation of helper tools

PAGAN can use MAFFT and RAxML/BppDistTree to infer the guidetree and Exonerate to speed up the alignments; it can also BppAncestor to compute maximum likelihood ancestral character states . Those to work, you need to have MAFFT, RAxML, Exonerate and BppSuite installed and either on the execution path or in the same folder with the PAGAN executable (PAGAN will use the latter if available.).  
 

Note that the PAGAN binary packages come with all the helper tools and, if they work on your system, they do not need to be installed separately.

You can download the programs from [http://www.ebi.ac.uk/~guy/exonerate/](http://www.ebi.ac.uk/~guy/exonerate/) (Exonerate), [http://mafft.cbrc.jp/alignment/software/source.html](http://mafft.cbrc.jp/alignment/software/source.html) (MAFFT) and [http://sco.h-its.org/exelixis/software.html](http://sco.h-its.org/exelixis/software.html) (RAxML), and follow the instructions for their installation. On many popular Linux distributions these programs are available in the software repository and can be installed pre-compiled. On Ubuntu, that is done with the following command:

```
sudo apt-get install exonerate
sudo apt-get install mafft
sudo apt-get install raxml
```

Because of dependencies to other programming libraries, the installation on MacOSX is easiest with a package management system such as [Homebrew](http://mxcl.github.com/homebrew/). On OSX 10.7, Exonerate and MAFFT can be installed with the following commands:

```
# uncheck below to install Homebrew
# sudo /usr/bin/ruby -e "$(curl -fsSL https://raw.github.com/gist/323731)"
sudo brew install pkg-config
sudo brew install exonerate
sudo brew install mafft
```

See [http://mxcl.github.com/homebrew/](http://mxcl.github.com/homebrew/) for further details.

BppSuite is provided as pre-compiled binaries for several platforms at [http://home.gna.org/bppsuite/](http://home.gna.org/bppsuite/) and is also available in repositories of some popular Linux distros. Unfortunately, the version of BppAncestor provided may contain a critical bug and may not work. See [here](https://github.com/ariloytynoja/prank-msa/blob/master/docs/prank_installation.md#installation-of-helper-tools) for the instructions to re-compile the BppSuite tools.

[back to top](#introduction)

* * *

### Advanced: installing boost from source code

If pre-packaged files are not available for your platform or you are not entitled to install packages and cannot persuade your system administrator to do that, the Boost libraries necessary for the compilation of the PAGAN source code can be installed using following commands:

```
# make a temporary directory
mkdir ~/tmp_boost
cd ~/tmp_boost

# get the source code
wget -O boost_1_55_0.zip http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.zip/download
# OR (on OSX)
# curl --location -o boost_1_55_0.zip http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.zip/download
unzip boost_1_55_0.zip
cd boost_1_55_0

# compile and install
sh ./bootstrap.sh --prefix=$PATH_TO_PAGAN_DIR/boost --with-libraries=program_options,regex,system,thread
./bjam threading=multi link=static --prefix=$PATH_TO_PAGAN_DIR/boost install
# on MacOSX, this may be preferred
# ./bjam macosx-version=10.7 threading=multi link=static --prefix=$PATH_TO_PAGAN_DIR/boost install

# clean up
cd ../..
rm -r ./tmp_boost
```

This compiles static versions of the boost libraries. The good thing is that those get statically linked to the executable and one does not need to worry about library paths. The bad thing is that the executable is rather big in size.