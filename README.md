# This image belongs to a larger project called Bioinformatics Docker Images Project (http://pegi3s.github.io/dockerfiles)
## (Please note that the original software licenses still apply)

This image facilitates the usage of rosettadock-screening, a tool to help find interactors of a given receptor, using RosettaDock. Because of legal restrictions, you need to build the pegi3s/rosettadock image with the latest tag youself (see instructions [here](http://bdip.i3s.up.pt/container/rosettadock_builder).

# Using the RosettaDock-screening image in Linux

In order to be able to use the pegi3s/rosettadock-screening image you must first build the pegi3s/rosettadock:latest image in your computer first, due to legal issues. The instructions are [here](http://bdip.i3s.up.pt/container/rosettadock_builder). Then, you should adapt and run the following command: `docker run -it -v /your/data/dir:/data -v /var/run/docker.sock:/var/run/docker.sock -e HOST_DATA_DIR="/your/data/dir" pegi3s/rosettadock-screening`

In this command, you should replace:
- `/your/data/dir` to point to the directory where the data and configuration files are located.

Ligands must be placed in a folder named ligands located on /your/data/dir
The receptor must be placed on /your/data/dir

A file named config with the following structure (lines marked with # may be ommited) must also be present on /your/data/dir

```
[Input_Files] 
# The PDB receptor file
pdb_receptor = 

[Regions]
# Region, Site, Weight
#For instance, (no line limit)
D185 = 1, 85, 0.80
D2190 = 2, 190, 0.51

[Variables_TF] 
# True and False centroid values for each region (no line limit)
T11 = 98
F11 = 19
T12 = 65
F12 = 10

[General_Constants] 
# Global variables for the pipeline (i stads for iteractions)
i_global_max = 100
i_local_max = 20
success_cycles = 10

[Reward] 
ratio = 1
ratio_increment = 0.1

[CST_Weight] #default is 0.1
cst_weight = 0.05
```
and you also need in /your/data/dir a constraints.cst file that looks like:

```
SiteConstraint CA 85 B SIGMOID 15.0 0.1
SiteConstraint CA 190 B SIGMOID 15.0 0.1
```
If a file named control_results.xlsx (optional) is present on /your/data/dir with data on the same ligands using a different method (Haddock, for instance), it will be used to produce a file that shows the results obtained with pegi3s/rosettadock-screening and the other method, for the same ligands, side by side. The structure of the .xlsx file must be for two regions "Ligand name", "Region 1", "Region 2", "Distance to positives", "Distance to negatives" 


