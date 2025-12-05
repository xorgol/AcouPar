# AcouPar 
AcouPar stands for "Acoustical Parameters" and it is the Command Line version of the xfm plugin "[ISO3382 Acoustical Parameters](https://www.angelofarina.it/Aurora/download/Aurora45-Alpha/acoustic45.xfm)" of the [Aurora plugins suite](http://www.angelofarina.it/aurora/download/Aurora45-Alpha/).

The plugin is meant for computing acoustical parameters according to ISO3382-1.
This command line utility produces exactly the same results of the plugin.

AcouPar processes a single or a couple (stereo file)  of Impulse Responses in a wavefile providing different approaches:
 1. One for single or two independent omnidirectional microphones impulse responses
 2. One for ne for WY-Ambix impulse responses
 3. One for WY-Ambix impulse responses (omni & figure-of-eight) pressure-velocity (pu)
 4. One for binaural impulse responses

Other approaches are available but not yet deeply testes such: 

 - One for p-p sound intensity probe
 - One for pressure-velocity with old FuMa WY

The command line utility is multiplatform and it can be run both on Windows-based machines and MacOs.
A pre-compiled version is already avaliable of it can be easily compiled from scratch.

### Output
The elaboration generate a txt file containing the acoustical parameters values.
The file is named after the stereo mode that has been chosen. For the 2 (single) omni elaboration the file will be named "AcouPar_omni.txt"; for the binaural instead will be named AcouPar_BIN.txt.


The inside of the file is arranged by parameters and per frequency in a tab-separated manner. In this way the data can be easily imported in a spreadsheet or other software for further elaboration.

**N.B.: decimal point is always used despite the region settings of your machine.**

```
Acoupar.exe - ISO3382 Acoustical Parameter File	
Filename	Par1_freq1  Par1_freq2  ..  Par2_freq1  ..
file.wav	value       value       ..  value       ..
```

Note also that if a file with the same name already exists, the output will be appended to the existing file.

## How-to Use
To use the Command Line tool a mono or stereo wavefile containing the Impulse Responses is required, and it shall be placed in the same directory of the binary. 

The stereo option must be specified and it is one of the following:
 - `--omni` for omnidirectional microphone(s)
 - `--wy` for soundfield analysis of WY-ambix
 - `--pu` for pressure-velocity analysis WY-ambix 
 - `--bin` for binaural analysis 
 
Then, the following code can be run on Windows systems:

```
.\AcouPar.exe --pu filename_WY.wav
```

While on MacOs:
```
./AcouPar --pu filename_WY.wav
```

A relative path to the wavefile can also be specified. Remember that in case that the filename or its path contains spaces, it is necessary to bracket it with double quotes like:
```
.\AcouPar_pu.exe "C:\folder with spaces\file name with spaces WY.wav"
```

## How-to Build
The CMake GUI can be employed to generate the necessary Visual Studio or XCode project data for building within the IDE. 

A quicker way requires to run the following commands from the inside of the repository folder.

```
mkdir build
cd build
cmake ..
```

then, for Windows systems run:
```
cmake --build .\
```

while for MacOs:
```
make
```

## Note

Please note that currently we are not supprting yet the two other possible formats for stereo
impulse responses: Soundfield WY (the old FuMa format) and p-p sound intensity probe.
These will possibly be added in the future.

## Copyright
(C) Angelo Farina, 15 November 2020