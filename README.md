---
# garfieldpp_nitpc
## Garfield++ toolk kit for negative ion simulation

---
# README Change Log
2021/02/16 T.Shimada: First commit
2021/06 H.Ishiura: Small modification. 
---
## Main source code author 
Takuya Shimada (Kobe Univ.)


## Maintaner 
Hirohisa Ishiura  (Kobe Univ.)

---
## garfieldpp_nitpc Overview
This code contains new class "AvalancheNIMicroscopic" for negative ion gas MPGD simulation in Garfield++.
Source/AvalancheNIMicroscopic.cc and Include/AvalancheNIMicroscopic.hh are added for negative ion tracking and MPGD simulation.




## Usage

# .bashrc
```
export GARFIELD_HOME=/pash/to/garfieldpp_nitpc
export HEED_DATABASE=$GARFIELD_HOME/Heed/heed++/databece
```

```
$ cd $GARFIELD_HOME
$ make
```
---
## Example
```
$ cd $GARFIELD_HOME/Example/NITPC/src
$ vim nitpc.cxx
```

nitpc.cxx
```
// nitpc.cxx
gas->LoadIonMobility("/path/to/Data/IonMobility_SF6-_SF6.txt"); // Mobility data path, user should set this file
...
std::string data_dir = "/path/to/gem/"; // Elmer output directory
std::string header = data_dir + "/gemcell/mesh.header";
std::string element = data_dir + "/gemcell/mesh.elements";
std::string node = data_dir + "/gemcell/mesh.nodes";
std::string eps = data_dir + "/gemcell/dielectrics.dat";
std::string volt = data_dir + "/gemcell/gemcell.result";
...
```

and,
```
$ make
$ cd ..
$ ./nitpc
....
create track.png
....
$ display track.png 
```
track.png shows the visualization of this simulation.

Time lapse animation
![avalanche_negativeion_gem](https://user-images.githubusercontent.com/52315643/108066478-bcd4af00-70a2-11eb-880f-beb927d8eee5.gif)

