# Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory

Python library for calculating Floquet vectors and simulating dynamics for fast entrainment by strong periodic input.

The codes in this repository implement the methodologies presented in 
S. Takata, Y. Kato, H. Nakao, Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory, Chaos: An Interdisciplinary Journal of Nonlinear Science, 31(9), 093124 (2021).

Please cite the paper if you find the codes are useful. 


## Integrated Development Environment(IDE) 

We use the Spyder IDE [https://www.spyder-ide.org/] to run the codes.

Note that the standard cell (section) separator (#%%) may not work in other IDEs.


## Setup 

Install requirements  

pip3 install -r requirements.txt


## Running the code

A brief description of the code relationships can be found in "code_connection.pdf".

### Stuart-Landau model

The main code is "SL_simple.py".   

#### SL_simple :  
This code simulates a Stuart-Landau model subjected to the optimal input proposed by Zlotnik (2013). The code in the first section calculates the Lagrange multipliers for the input, simulates the dynamics with the input, and stores the simulation results in arrays. The code in the second section draws the figure and saves the simulation results. The figure plotted by the code corresponds to Fig.1 in the paper.

### van der Pol model

The main codes are "VAN_Floquet", "VAN_simple", "VAN_feedback", "VAN_penalty", "VAN_tangent", "VAN_plot", and "calculate_arnold_tongue".   

#### VAN_Floquet :   
This code calculates the limit cycle, frequency, Floquet vectors, etc., of a van der Pol model. The code in the first section calculates the above quantities. The code in the second section draws the figure and saves the simulation results. The figure plotted by the codes corresponds to Fig.2 in the paper.
  
#### VAN_simple : VAN_Floquet needs to be run beforehand.
This code simulates a van de Pol model subjected to the optimal inputs proposed by Zlotnik (2013). The code in the first section calculates the Lagrange multipliers for the input. The code in the second section simulates the dynamics with the input and stores the simulation results in arrays. The code in the third section calculates the phase coupling function and saves the simulation results.
  
#### VAN_feedback : VAN_Floquet needs to be run beforehand.
This code simulates a van der Pol model subjected to the input of the amplitude-feedback. The structure of the sections is the same as "VAN_simple".
  
#### VAN_penalty : VAN_Floquet needs to be run beforehand. 
The code simulates a van der Pol model subjected to the input of the amplitude-penalty method. The structure of the sections is the same as "VAN_simple".
  
#### VAN_tangent : VAN_Floquet needs to be run beforehand. 
The code simulates a van der Pol model subjected to the input of the tangent-only method. The structure of the sections is the same as "VAN_simple". 
  
#### VAN_plot : It is necessary to save simulation data to draw figures in advance.  
This code draws figures for the simulation results of van der Pol models. In the first section, we set the matplotlib parameters. In the second and third sections, we draw and save the simulation results with the amplitude-feedback and amplitude-penalty methods, respectively. The figure plotted in the section corresponds to Fig. 3 in the paper. In the fourth section, we draw and save the simulation results to compare the two proposed methods and the method proposed by Zlotnik (2013). The figure plotted in the section corresponds to Fig. 4 in the paper. In the fifth section, we draw and save the simulation results of the tangent only method. The figure plotted in the section corresponds to Fig. 8 in the paper.
  
#### calculate_arnold_tongue : VAN_Floquet needs to be run beforehand.
This code simulates three methods, draws the figure of the Arnold Tongue, and saves the simulation results. In the first section, we set up the three functions needed in the calculation. In the second, third, and fourth sections, we simulate the dynamics, draw the figure of the Arnold tongue and save the results with the optimal input proposed by Zlotnik (2013), the amplitude feedback method, and the amplitude penalty method, respectively. The figure plotted in the section corresponds to Fig. 5 in the paper. In the fifth section, we draw and save the color bars.

### Willamowski-Rossler

The main codes are "WR_Floquet", "WR_simple", "WR_feedback", "WR_penalty", and "WR_plot". The structure of these codes is the same as the corresponding codes of "VAN_...". The figures plotted in the code are corresponding to Fig. 6 and 7 in the paper.

### FitzHugh-Nagumo

The main codes are "FHN_Floquet", "FHN_simple", "FHN_feedback", and "FHN_penalty". The structure of these codes is the same as the corresponding codes of "VAN_...". Note that the results for the FitzHugh-Nagumo model are not presented in the paper. 


## Questions 

Please feel free to contact me (tktsho72@gmail.com) if you have any questions.
