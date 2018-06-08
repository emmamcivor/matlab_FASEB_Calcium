# Instructions for running code

This code was used to create plots of the simulated calcium profiles in the poster 'Assessing the impact of Orai1 channel cross-linking on store operated ER refilling: Insights from a mathematical model' presented at the FASEB Calcium and Cell Function conference in Tahoe City, June 2018.

Below are the instructions for running the code, everything you need is within this folder. 
All instructions are run in a MATLAB session, and assume that the scripts are available in the working directory.
There are a selection of sample simulations, these are the simulations from the poster but are completely independent of each other so can be run individually to create the plots used to make the figures on the poster. 

The code used to create the figures in the poster used a fine spatial grid (dx=2nm, dy=2nm), as in file

	soce_js_HEK_config.m

This takes about 30 minutes to run on one CPU core of a Dell R630 server with 2 x Intel(R) Xeon(R) E5-2660 v4 processors and 768Gb RAM. 

I have included code with a coarse spatial grid (dx=5nm, dy=5nm), in file

	soce_js_HEK_config_coarse.m

This takes about 5 minutes to run on a laptop.

The code is a prototype and as part of our subsequent work we intend to develop more user friendly software, for now if you have any problems or would like help please contact us. 
A word about memory usage: the code creates large matrices and the simulations are saved as four dimensional matrices, approximately 2GB in size; these are all held in memory. 

## Example simulations

These are six separate independent example simulations. Each example can be run individually and will plot the relevant calcium profiles used to create the figures in the paper. I have provided instructions for the coarse simulations, file 

	soce_js_HEK_config_coarse.m

but the same instructions should be applied for the fine simulations using file

	soce_js_HEK_config_coarse.m

If you want to save the output from the simulations edit the file

	soce_js_HEK_config_coarse.m

or 

	soce_js_HEK_config.m

and change the following parameter as shown.

	save_data=1;

This will save three dimensional matrices totalling approximately 1Gb per simulation.

### Cross-linked Orai channels surrounded by a ring of SERCA pumps 30nm away

Edit the file

	soce_js_HEK_config_coarse.m

by changing the following parameter values as shown and save your changes.

	situation=1;

This ensures the Orai channels are cross-linked and the distance between the Orai channels and SERCA pumps is 30nm. 

Execute the following command.

    run soce_js_HEK_config_coarse.m

This will simulate SOCE and the calcium diffusion and ER refilling occurring in 1ms and then plots the calcium profiles in Figs. 4 and 6.

### Non cross-linked Orai channels surrounded by a ring of SERCA pumps 30nm away

Edit the file

    soce_js_HEK_config_coarse.m

by changing the following parameter values as shown and save your changes.

    situation=2;

This ensures the Orai channels are not cross-linked and the distance between the Orai channels and SERCA pumps is 30nm.

Execute the following command.

    run soce_js_HEK_config_coarse.m

This will simulate SOCE and the calcium diffusion and ER refilling occurring in 1ms and then plots the calcium profiles in Figs. 5 and 6.

### Cross-linked Orai channels with SERCA pumps placed along the periphery of the ER-PM junction

Edit the file

    soce_js_HEK_config_coarse.m

by changing the following parameter values as shown and save your changes.

    situation=3;

This ensures the Orai channels are cross-linked and the SERCA pumps are along the periphery of the junction.

Execute the following command.

    run soce_js_HEK_config_coarse.m

This will simulate SOCE and the calcium diffusion and ER refilling occurring in 1ms and then plots the calcium profiles in Fig. 7.

### Non cross-linked Orai channels with SERCA pumps placed along the periphery of the ER-PM junction

Edit the file

    soce_js_HEK_config_coarse.m

by changing the following parameter values as shown and save your changes.

    situation=4;

This ensures the Orai channels are not cross-linked and the SERCA pumps are along the periphery of the junction.

Execute the following command.

    run soce_js_HEK_config_coarse.m

This will simulate SOCE and the calcium diffusion and ER refilling occurring in 1ms and then plots the calcium profiles in Fig. 7.

### Cross-linked Orai channels with SERCA pumps placed along the periphery of the ER-PM junction and 12 extra SERCA pumps in ER-PM junction

Edit the file

    soce_js_HEK_config_coarse.m

by changing the following parameter values as shown and save your changes.

    situation=5;

This ensures the Orai channels are cross-linked, the SERCA pumps are along the periphery of the junction and extra SERCA pumps are included.

Execute the following command.

    run soce_js_HEK_config_coarse.m

This will simulate SOCE and the calcium diffusion and ER refilling occurring in 1ms and then plots the calcium profiles in Fig. 8.

### Non cross-linked Orai channels with SERCA pumps placed along the periphery of the ER-PM junction and 8 extra SERCA pumps in ER-PM junction

Edit the file

    soce_js_HEK_config_coarse.m

by changing the following parameter values as shown and save your changes.

    situation=6;

This ensures the Orai channels are not cross-linked, the SERCA pumps are along the periphery of the junction and extra SERCA pumps are included.

Execute the following command.

    run soce_js_HEK_config_coarse.m

This will simulate SOCE and the calcium diffusion and ER refilling occurring in 1ms and then plots the calcium profiles in Fig. 7.


