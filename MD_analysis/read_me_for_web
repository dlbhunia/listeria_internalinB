This file contains example data that was analyzed from the full trajectory described in (Liu et al., 2023).
The data here only represents the first few ns of a 200 ns run. The goal of this example is to allow users to run the code that was generated to analyze the files on a much smaller and user-friendly set of data. The data will not look the same as the manuscript as a result. 
This folder should have the following contents:
1)	Read_me.docx
  a.	This document containing instructions on how to run the code
2)	code
  a.	Chain_distances.py
    i.	File used to analyze the trajectory in UCSF Chimera
  b.	Distance_Analysis.py
    i.	File used to analyze the outputs from the trajectory analysis
3)	MD_files
  a.	run_input.tpr
    i.	topology and run parameters for the MD simulation
  b.	trajectory.xtc
    i.	trajectory file for the MD simulation
4)	output
  a.	example.txt
    i.	example output file generated from Chain_distances.py file. File is analyzed using the Distance_Analaysis.py file
  b.	norm_example.png
    i.	example image that will be generated from the code
    Formatting notes
    >>> is the symbol for a line of code to run. This could be a command line prompt or Python code.
    > Button to press in software. For example, to turn your computer on it would be written:
    >power

To use this file:
  1)	Install required software
    a.	Chimera 
      i.	UCSF Chimera, version 1.15 was used in this project. This version is no longer supported by UCSF Chimera, which may lead to the analysis breaking if run on newer versions. This is especially true of analyzing files using Chimera X. 
      ii.	https://www.cgl.ucsf.edu/chimera/
    b.	Python
    i.	Install Python. This was written in Python 3.8.12
    ii.	https://www.python.org/
  c.	Anaconda
    i.	Install Anaconda. This was written in Anaconda 22.9.0.
    ii.	https://anaconda.org/
  d.	Spyder
    i.	Install Spyder through the Anaconda platform. This was written using Spyder 5.2.1.
  e.	Libraries
    i.	It is recommended that you create a virtual environment in Anaconda before conducing any analysis. In the new environment, install pandas and matplotlib. 
    ii.	Pandas version 1.3.5
    iii.	Install pandas
        >>> conda install pandas
    iv.	Install matplotlib 3.5.1
        >>> conda install matplotlib
2)	Analyze molecular dynamics trajectory
  a.	Files required
    i.	run_input.tpr
    ii.	trajectory.xtc
    iii.	Chain_distances.py
  b.	Files generated
    i.	Example.txt
  c.	Place the required files in a location on your computer where you can access them. 
    i.	It may be convenient to make a .txt file with frequently used file paths
    ii.	In all cases below where PATH\TO\file_name.type is used, it indicates you should replace that with the file path on your own device. 
  d.	Change file path in the Chain_distance.py file
    i.	Change the 14th line of the Chain_distance.py file to a folder where you want your MD simulation results to be saved. This can be done in Spyder, or any other IDE that you prefer. 
    ii.	If you are analyzing the MD files contained above, this program will give you the LAP residues that are within 5Å of internalin. If converting to a different trajectory, the code can be modified by changing the residues that are selected on line 26 and 36 and renaming selections to more intuitive names based on your project in lines 28 and 32. 
  e.	Open the simulation in UCSF Chimera. 
    i.	Open UCSF Chimera
    ii.	Open the simulation
    1.	>Tools > MD / Ensemble Analysis > MD Movie
      iii.	In the window that pops up use the following inputs:
          1.	Trajectory format = GROMACS
          2.	Run input (.tpr) = PATH\TO\run_input.tpr
          3.	Trajectory (.trr or .xtc) = PATH\TO\trajectory.xtc
          4.	Use frames = first through =last
      iv.	This should open the simulation:
 
    1.	Your Chimera display may look different depending on what you have in the quick access tool bar on the left-hand side. The snip also shows the movie player in the window, which can be moved to any convenient location.
      v.	Play the full movie. The file that you have opened is just the first few ns of the simulation from (Liu et al., 2023) (26 frames ~1 ns). This loads all of the frames you have opened. 
      f.	In the command line type “open”
        i.	Use the browser window to select the Chain_distances.py file
        ii.	When you click open in the browser on this file, it will start to run on the program. This may take a few minutes depending on the computer and the size of your simulation. If will appear that Chimera is frozen while the program runs. It is still running though, so be patient.
      g.	If successful, a new file will be generated in the path your outputs names “example.txt” or another name if you chose to change the file name. This is the file that is analyzed with the next Python code in this example.
3)	Analyze amino acids results from Chimera Distance Analysis
  a.	Files required
    i.	Distance_Analysis.py
      1.	This file is used to analyze individual amino acid counts from the trajectory and group them based on charge and hydrophobicity, and it does it for the interaction B and C chains separately. It was written to count each amino acid type individually instead of in groupings by default, since many combinations of amino acid groupings were considered during the original data exploration. This could be updated and made significantly better, but its current format is the result of the data exploration process. 
  b.	Files generated
    i.	norm_exmple.png
  c.	Open Anaconda
  d.	Open Spyder
  e.	From Spyder, open Distance_Analysis.py
  f.	In this file, each cell is indicated by #In[CELL NAME]. A summary of cell functions follows. These cells can be run individually, or the whole program can be run to visualize different results from the analysis. What follows is a brief description of each cell. 
  g.	#In[IMPORT LIBRARIES] 
    i.	Read in the dependencies of this project including Pandas and Matplotlib
  h.	#In[READ AND CLEAN THE DATA]
    i.	This section reads in the file generated from the Chain_distances.py analysis file
    ii.	It strips unwanted punctuation that is an artifact of the program that writes the data and assigns a time value to each step. 
    iii.	The data format is now one row for each frame in the MD simulation and a series of columns that show which amino acids were within 5Å of internalin in the columns.
  i.	# In[MAKING A DF OF A PARTICULAR CHAIN IF NEEDED]
    i.	This looks at the interacting B and C chains separately. It may not be necessary for your project depending on the interactions you are interested in investigating. 
  j.	# In[MASKING RESIDUES TO PLOT SPECIFIC CATEGORIES]
    i.	This section finds the counts of interacting residues for each time, by type, for each frame in the MD simulation. 
    ii.	Amino acids are grouped into different categories including positive, negative, and hydrophobic (Kyte-Doolittle hydrophobicity > 0).
    iii.	Depending on the size of your trajectory, this may take some time to complete, as the process is conducted for each amino acid individually. It could be sped up by putting parts of this into functions and pre-grouping amino acid categories. 
  k.	# In[APPLY A ROLLING MEAN]
    i.	Applies a rolling mean to the different interaction counts calculated above. Since this examples data is fairly small, a rolling mean of 10 points was selected, which is different than what was originally published, where the rolling mean was 300 points. 
  l.	# In[GENERATE THE NORMALIZED VALUES]
    i.	This section normalizes that counts to their frequencies in the B and C chains. The chain length and amino acid / amino acid type counts are the sum of those amino acids from the B and C chain in the LAP model.
  m.	# In[PLOTTING INDIVIDUAL INTEACTIONS- COUNTS]
    i.	This section shows how the individual groupings of amino acids can be plotted. These can be updated to fit your project needs. 
  n.	# In[PLOTTING ALL INTERACTIONS- COUNTS]
      i.	This section plots all of the counts of amino acids on top of one another. It can be seen that there is not much distinction between interactions. 
  o.	# In[PLOTTING ALL INTERACTIONS- NORMALIZED]
    i.	This section plots the normalized data. This shows that there is a distinction between interacting amino acids when they are normalized to their frequencies within the chains.  
    ii.	This figure is also saved for an example. 

