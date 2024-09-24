# PharmaSuiteProject
This was my final project for an Intro to Python for Life Sciences Class (BME 160) at UCSC. This is a work in progress and ostensibly does not run since it requires such an intense effort on the part of anyone wishing to use it to configure it. However I have coerced it previously into producing results which can be seen in the test results PDF. I promise although it is a mess I have reason to be proud enough to upload it publically. The concept of this software is to rank drug candidates from .sdf files produced by several companies who dump massive lists of theoretical drug candidates for treatment of various diseases. The software reads and writes .sdf files and ranks candidates contained therein based off of their theoretical bioavailability according to some calculations derived from a few chemical properties. The real promise of the software is the docking component. I was shocked to see looking into softwares to model ligand binding that there were very limited options to automatically dock large groups of ligands. I expected this to be the easy part of the project and had many other functionalities planned but the meat of the code ended up being dedicated to creating what I believe to be the first program capable of automatically scoring ligand binding affinity to enzymes using a large group of ligands stored in a .sdf file. This involves the arduous process of converting the coordinate data stored in the .sdf files to the .pdb format and then converting that to the .pdbqt format and then finally of using the command functionality of the autodock vina software and some batch scripts to dock and store the docking ratings of all the candidates. The example application we presented was a list of candidates for tuberculosis drugs which bind to an enzyme humans also possess a version of. The theory was to rank candidates on their bioavailability and also on their potential to bind to this concerning off target site that might induce seizures or other side effects. We did manage to rank a list of some 100 odd candidates and get an A-, which I am in hindsight quite grateful for given that the professor had to take it on faith that the code worked based off visual inspection and the charts we produced.
