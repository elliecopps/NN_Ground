# SGViz_InfD
Visualization of Infinite Dimensional Spin Glass

## Getting it working
To run this you'll need to be running some version of Python 3. Additionally you'll need to install the PyQt library with 
```
pip3 install pyqt5
```
You'll also need to make sure you have the libraries at the top of the ```SGViz.py``` file. 

## Running the code 
 All you need to do is run
```python3 SGViz.py```

## The interface
The interface is pretty intuitive. Enter values in the text boxes and press the Enter button or press Enter on your keyboard to implement the change. 

* N - Number of spins, integer values only, gets slow past 16 or so

* Seed - seed of the Jij Matrix, integer values only

* Config - Configuration of the spins, integer values only from 0 to 2^N - 1 , 0 is all spins down 2^N - 1 is all spins up

* Bipartition<sup>*</sup> - choice of bipartition to split spins into two parts (equal parts if N is even), integer values only from 0 to total number of possible bipartitions 1/2 ( N choose N/2 )

* The Clear button - clears everything off the visualization except the visualization, I used this when taking screenshots. Press Clear again to bring everything back.

Additionally, in the top right, the current configuration is shown in terms of 1s and 0s. 0 is spin down, 1 is spin up. Moreover, the configurations that have the lowest energy are shown in the GS Configurations array. The current energy of the configuration is shown, as well as the total number of satisfied and unsatisfied bonds. 

<sup>*</sup> This functionality initially turned off because it's not always useful. See the next section to turn it on

## Changing things

I haven't implemented command line arguments. The two functionalities that I thought would be useful to turn off and on (but haven't created buttons for them yet) are toggling whether to show bipartitions, and whether to show the ground state array and current energy. The former is because the bipartition visualization is annoying if you don't want to see it. The latter is because it's computationally expensive to calculate the ground state array every time. If you open ```SGViz.py``` you'll see the following code near the top

```
#Set these to True to see bipartitions or Energy info, respectively
######################################################################
self.show_bps = False
self.show_Ens = True
######################################################################
```

Change the first boolean to True to see the bipartitions, and the second boolean to False to stop showing the ground state array and current energy.
