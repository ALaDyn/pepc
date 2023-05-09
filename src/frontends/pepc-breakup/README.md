# Usage of pepc-breakup frontend

After the `make pepc-breakup` command ran succesfully, a folder named `bin` should be created with the compiled executable within. 

Copy the files from `src/frontends/pepc-breakup/bin_files/` into the `bin` folder. 
Create a copy of `params_template` and `job_template` respectively, adapt the copied files to your specific requirements.

In the case of tokamak simulation with a set of multiple poloidal coils, create your own text file which describes the 
dimension and the current flowing through each coil. This determines the poloidal magnetic field structure within the tokamak domain.

The relative path to the created coil description file should be given to the `coil_data_file` variable to enable correct computation
of the poloidal magnetic field structure.
