# Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory

Python library for calculating floquet vectors and simulating dynamics.   

The codes in this repository implement the methodologies described in   
Takata, Shohei, Yuzuru Kato, and Hiroya Nakao. "Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory." arXiv preprint arXiv:2104.09944 (2021).  

Please cite this paper when you use the codes.  

## Integrated Development Environment(IDE) 

The code is intended to be run on spyder[https://www.spyder-ide.org/].

Please note that the section marks (#%%) will not work when run in other IDEs. 


## Setup 


Install requirements  

pip3 install -r requirements.txt


## Running the code

A brief description of the code relationships can be found in "code_connection.pdf".

### Stuart-Landau

The main code is "SL_simple.py".   

#### SL_simple :  
This code simulates the optimal inputs by Zlotnik(2013)for Stuart-Landau model.
The first section calculates the Lagrange multipliers for the input, simulates the dynamics with the input, and stores the data in arrays.
In the second section, we draw and save the simulation results. The figure corresponds to fig1 in the paper*.

### van der Pol

The main codes are "VAN_Floquet", "VAN_simple", "VAN_feedback", "VAN_penalty", "VAN_tangent", "VAN_plot", and "calculate_arnold_tongue".   

#### VAN_Floquet :   
This code calculate limit cycle, frequency, Floquet vectors, etc. for van der Pol model.  
The calculation is done in section 1, the drawing of the figure and saving of the figure and data is done in section 2. The figure corresponds to fig2 in the paper*.
  
#### VAN_simple : It is necessary to run VAN_Floquet beforehand.  
The code simulates optimal inputs by Zlotnik(2013) for van der Pol model.   
The first section computes the Lagrange multipliers for the input.  
In the second section, we simulate the dynamics with the input and keep the simulation results in arrays.  
In the third section, we calculate the phase coupling function and store the simulation data.
  
#### VAN_feedback : It is necessary to run VAN_Floquet beforehand.  
This code simulates amplitude-feedback method for van der Pol model.   
The structure of section is the same as "VAN_simple".  
  
#### VAN_penalty : It is necessary to run VAN_Floquet beforehand. 
The code simulates amplitude-penalty method for van der Pol model.   
The structure of section is the same as "VAN_simple". 
  
#### VAN_tangent : It is necessary to run VAN_Floquet beforehand. 
The code simulates tangent-only method for van der Pol model.   
The structure of section is the same as "VAN_simple"  
  
#### VAN_plot : It is necessary to save simulation data to draw figures in advance.  
This code draws simulation results of van der Pol model.   
In section 1, we set the matplotlib parameters.  
Section 2 draws and saves the simulation results of the amplitude feedback method. The figure corresponds to fig3 in the paper*.  
Section 3 done for the amplitude-penalty method. The figure corresponds to fig3 in the paper*.  
Section 4 done for the amplitude-feedback and amplitude-penalty methods. The figure corresponds to fig4 in the paper*.  
In section 5, we draw and save the simulation results for periodic input in the tangential direction. The figure corresponds to fig8 in the paper*.  
  
#### calculate_arnold_tongue : It is necessary to run VAN_Floquet beforehand.  
This code simulates three methods, draws the Arnold Tongue, and saves the figure and data.  
Section 1 sets up the three functions for the next sections.  
In section 2, we simulate the Zlotnik optimal input (2013), draw the Arnold tongue, and save the figure and data. The figure corresponds to fig5 in the paper*.  
In section 3, we do the same with the amplitude feedback method.  
In section 4, we do the same with the amplitude penalty method.  
Section 5 is about drawing and saving the color bars.  

### Willamowski-Rossler

The main codes are "WR_Floquet", "WR_simple", "WR_feedback", "WR_penalty", and "WR_plot".   
The structure of these codes are the same as "VAN_...".  
The figure is corresponding to fig6 and fig7 in the paper*.  


### FitzHugh-Nagumo

The main codes are "FHN_Floquet", "FHN_simple", "FHN_feedback", and "FHN_penalty".   
The structure of these code are the same as the above.  
This model is not attached in the paper*.   
  
  
*Takata, Shohei, Yuzuru Kato, and Hiroya Nakao. "Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory." arXiv preprint arXiv:2104.09944 (2021).


## Questions 

If you have any questions, don't hesitate to email takata.s.ae@m.titech.ac.jp
