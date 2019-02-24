# Aggregate
python post processing tool based on vtk and pandas to derive aerodynamic coefficients from triangular surfaces input(s) 

## Table of contents
* [TODO](#todo)
* [Usage](#usage)
   * [Integrate](#integrate)
   * [Concatenate files](#concatenate-files)
   * [Visualize data](#visualize-data)
   * [Complete exemple](#complete-exemple)
* [Dependancies](#dependancies)


## TODO
 - Implement and test "reverse_aero_convention". Only alpha=beta=0 has been tested
 - Add test 
 - Finish readme and documentation (transformation, data extraction using a serie of point) 

## Usage

```
usage: core.py [-h] -i LIST_INPUT [LIST_INPUT ...] [-f FORCEFORMAT]
               [-o OUTPUTNAME]
               [-k LIST_VARIABLETOKEEPFORWRITING [LIST_VARIABLETOKEEPFORWRITING ...]]
               [-Cp CP] [-Cfx CFX] [-Cfy CFY] [-Cfz CFZ] [-q Q] [-qd QD]
               [-alpha ALPHA_DEG] [-beta BETA_DEG] [-cog xcog ycog zcog]
               [-Sref SREF] [-Lref LREF] [-reverse_normal]
               [-reverse_aero_convention] [-autoOrient] [-c] [-t tx ty tz]
               [-r r11 r12 r13 r21 r22 r23 r31 r32 r33] [-v] [-nw]
               [-glvar VARIABLETODISPLAY] [-gl] [-e EXTRACTFROMCVS]
```


### Integrate
Define input data files:
 -  "-i"

Define variable to be used for integration :
 -  "-Cp"
 -  "-Cfx" 
 -  "-Cfy"
 -  "-Cfz"
 -  "-q"  

Reference variables:
 - "-Sref" (defaut = 1)
 - "-Lref" (defaut = 1)
 - "-cog" : center of gravity for moment (defaut = 0. 0. 0.) 

Orientation of the inflow :
 - "-alpha" : angle of attack [deg] (defaut = 0)
 - "-beta" : angle of attack [deg] (defaut = 0)

Orientation of the inflow :
 - "-alpha" : angle of attack [deg] (defaut = 0)
 - "-beta" : angle of side-slip [deg] (defaut = 0)

```
python3 Aggregate/core.py -i data/xcore-F1000-0.vtk -Cp "Cp [-]"  -Cfx "Cfx [-]" -Cfy "Cfy [-]" -Cfz "Cfz [-]"  -q "Heat Flux: Net [W/m2]"  
```

```
...
write integration data in : output.csv
plot in : output.png
```


![integrate1](./images/integrate1.png)


### Concatenate files
```
python3 Aggregate/core.py -i  data/xcore-F1000-0.vtk data/xpanel*  -c 
```
```
...
concatenate the inputs...
writting: output-concat.vtk
```
### Visualize data
Read files, concatenate them and visualize the results in an interactive windows.

#### exemple without specified the variable to be displayed 
```
python3 Aggregate/core.py -i  data/xcore-F1000-0.vtk data/xpanel*  -c -gl 
```
#### exemple selectioning the variable to display (auto-scaling)  
```
python3 Aggregate/core.py -i  data/xcore-F1000-0.vtk data/xpanel*  -c -gl  -glvar "Cp [-]" 
```

The output.vtk file can also be viewed using Paravie <url>www.paraview.org</url>
#### output
```
Available data
                BC [-]
                Cf [-]
                Cfx [-]
                Cfy [-]
                Cfz [-]
   SELECTED ->  Cp [-]
                Cq: Net [-]
                Normals
                Pressure [Pa]

```
![exemple1](./images/exemple1.png)

  
### Complete exemple 
TODO

### Data extraction 
Using a serie of 3D points as input, recover the data that is the closest for each point
TODO

## Dependancies

```
pip3 install -r requirements.txt
```
