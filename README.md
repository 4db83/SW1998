## SW1998
Stock and Watson (1998) Replication files (see the [SW1998](SW1998.pdf) pdf file for details on the paper). 
This runs on Matlab version R2023a with the following toolboxes installed.

### MATLAB
- Econometrics Toolbox                   
- Optimization Toolbox                   
- Parallel Computing Toolbox             
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox   

PS. The Parallel Computing & Symbolic Math Toolbox are not needed. Most likely it is only the 
Optimization Toolbox that is required. 

Zip download or clone the SW1998 directory to your local machine. Then, simply run the '*SW1998_MUE_replication.m*' by pressing F5 (or the run button), and it should work.
Let me know if there any problems.

### Ghostscript (for PDF figures)
- Install Ghostscript for Windows from https://ghostscript.com/download/.  
- Edit `localFunctions/print2pdf.m` so that `GSversion` and the `GS` path point to your installed `gswin64c.exe` (or `gswin32c.exe`).  
- Without Ghostscript the script still saves EPS files, but PDF copies will not be produced.


db, sthlm, 28.10.2025.
