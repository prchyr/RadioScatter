 <img src="doc/RSLogo.png" alt="RadioScatter"> 

This is the RadioScatter code. It is used to simulate radar reflections from ionization deposits. It is meant to run independently of, or within, a monte-carlo program (such as GEANT4). More documentation can be found in the included user manual.

### Prerequisites

these are the general prerequisites, please see below for system-specific notes


ROOT 6.08 or higher (https://root.cern.ch/downloading-root)


the module slac_rf is g4 source code that will produce histograms of voltages received in a scattering experiment. to run inside of GEANT4, you'll need to have that installed

GEANT4 (https://geant4.web.cern.ch/geant4/support/download.shtml)
The GEANT4 software is fully independent of RadioScatter, and its license is available here: (https://geant4.web.cern.ch/license/LICENSE.html). Per that license: "This product includes software developed by Members of the Geant4 Collaboration ( http://cern.ch/geant4)."

NOTE: When installing GEANT4, be sure to install with the install data option (-DGEANT4_INSTALL_DATA=ON) to be sure that things work.

you can use any monte-carlo program you want with RadioScatter, but this package includes some GEANT4 programs, built on their examples, that harness the full power of GEANT to make realistic rf scattering signals. if you want to use those, i'm gonna assume that you know how to use GEANT4 and ROOT.

#### system-specific notes:

on Ubunutu 20.04 LTS, GSL is not installed by default, so ROOT will not build MathMore, which it needs. So, before installing ROOT, install the libgsl-dev library and be sure that ROOT is built with MathMore support (should be the default behavior). 


************
## INSTALL
************

**********************
### radio scatter default installation
**********************
first be sure to have ROOT  installed. there's lots of documentation on building those available elsewhere. and if you want to use it with GEANT, you'll need to have that installed too. then:

```
cd /your/favorite/source/dir  (like /usr/local or something)
git clone https://github.com/prchyr/RadioScatter.git
cd RadioScatter
./install.sh
```

this will install the header files and shared libraries (for root analysis) to the right places. it will put the libs inside of /usr/local/lib and the header files in /usr/local/include. It will make a build directory (inside the RadioScatter directory) called build, where the build products are placed before they are installed in the proper location.

*********************************
### slightly fancier installation (or if you can't write to /usr/local) 
*********************************

say you want your libraries and headers installed in some other place than /usr/local, you can set the RS_INSTAL_DIR variable in your bashrc. then libraries will install to RS_INSTALL_DIR/lib and headers to RS_INSTALL_DIR/include

```
cd /your/favorite/source/dir
git clone https://github.com/prchyr/RadioScatter.git
export RS_INSTALL_DIR=/your/preferred/install/dir
cd RadioScatter
./install.sh 
```

tested on:

Ubuntu 16.04/gcc 5.4

Red Hat Enterprise server 6.9/gcc 4.8.5

Ubuntu 18.04/gcc 7.4



*****************************************************************************************
## using
*****************************************************************************************



the theory of operation is outlined in this paper:

http://arxiv.org/abs/1710.02883

Please see the User Manual for detailed documentation of how RadioScatter works, and the function definitions. 


how it works, breifly:

The RadioScatter module calculates rf scattering from an energy deposit in 4 space. It uses the position, energy deposited, region over which the energy was deposited, and the ionization energy of the material to get a number of free charges over which to perform the scatter, and a number density to get material effects (plasma screening, eave damping, etc.) it is most useful to have the deposits be from some monte-carlo generation, like GEANT4. you can probably fit it into any MC program with relative ease, or generate your own plasma cloud to scatter radio from. 

**********************************************
### some examples:
********************************************

the examples directory has some use cases that don't involve geant4. to compile, be sure to link to the radio scatter library, i e -lRadioScatter, and also the root libraries, like so:

```g++ -std=c++0x scatterFromCascade.cc -o scatterFromCascade `root-config --cflags --glibs --libs` -lRadioScatter ```

in the examples, there is a use case where a shower has been prepared using GEANT4 and then this is used in a standalone program to scatter radio, as an example of runnning radioscatter outside of a monte carlo program. various options are shown and explained in the example

Once compiled, you would then run the above example via

```./scatterFromCascade <plasma lifetime [ns]> <frequency [GHz]> <tx power [W]>```

*************
### GEANT4 usage
*************

there is a GEANT4 module called slac_rf, so to install it, do

```
cd /your/GEANT/app/dir
mkdir slac_rf_build
cd slac_rf_build
cmake /path/to/RadioScatter/example/GEANT4/slac_rf
make -b -j4
```

which will make a build directory somewhere you want it to live, then install the GEANT4 simulation there linked to the RadioScatter libraries, once you've installed them as above.

how it works within GEANT4:

the Radioscatter.hh file is a header file that needs to be included in the GEANT4 source. there should be several different ways to use it, but one way is to globally declare the RadioScatter object in the main simulation .cc file of the simulation. it is then best to pass the constructed radio scatter object to the action initialization and the run action and stepping action files of the simulation. this will allow access to different methods that are useful at different points in the sumulation lifetime. for example, the data is filled during the stepping action phase, but the ROOT histograms are cosed at the end of the run action, through a built in GEANT4 method called EndOfRunAction. 

there are also several commands available through the RSmessenger.hh/cc files which allow for manipulation of the radioscatter configuration through the GEANT4 macro files. things such as

/RS/setTxPos 5 0 3 m

and

/RS/setTxPower 100

/RS/setPlasmaLifetime 1

which can be used without the need to re-compile the GEANT4 source, making running simulations much more simple.

ok so to actually run the program,

cd /path/to/geant4/program/install/dir/  (default is RadioScatter/slac_build unless you provided an argument to install.sh)

./slac_rf -m run1.mac -f "/path/to/file/and/filename.root"

to set the root file name, (-f) which is optional. each run is stored in an entry in a tree. there is a default macro

run1.mac:
	-this will set the tx position and the energies and such, and then will simulate the received signal in 3 antennas


### Analyzing the root files:

To analyze the root files, you need to be sure to add the newly created libRadioScatter to your rootlogon.C. Add the following:

gROOT->ProcessLine(".L /usr/local/lib/libRadioScatter.so");

gROOT->ProcessLine("#include <RadioScatter.hh>");

which should be enough to load in events and plot them, etc. 

## Documentation

a user manual in pdf form can be found at doc/userManual.pdf, and it will attempt to be current.  you can generate one yourself via doxygen (if you have it installed) by running 

```


doxygen doc/Doxyfile
```

from the RadioScatter directory. 

enjoy. 
