# Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory

Python library for calculating floquet vectors and simulating dynamics. 

The code in this repository was prepared to implement the methodologies described in 
Takata, Shohei, Yuzuru Kato, and Hiroya Nakao. "Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory." arXiv preprint arXiv:2104.09944 (2021).


#############################################
## Integrated Development Environment(IDE) ##
#############################################

We assume that the code runs on spyder[https://www.spyder-ide.org/].  

If you run the code on other IDE, please note that the section mark (#%%) doesn't work. 


###########
## Setup ##
###########

Install requirements

pip3 install -r requirements.txt


######################
## Running the code ##
######################

The brief description of the code relationships is in "code_connection.pdf". 

# Stuart-Landau

The main code is "SL_simple.py". 

SL_simple :
The code simulates Zlotnik optimal inputs(2013) for Stuart-Landau model. 
The first section is calculating lagrange multipliers for inputs, simulating dynamics with the input, and saving data.
The second section is drawing and saving simulation results. The figure is corresponding to fig1 in the article*.


# Van der Pol

The main code is "VAN_Floquet", "VAN_simple", "VAN_feedback", "VAN_penalty", "VAN_tangent", "VAN_plot", and "calculate_arnold_tongue". 

VAN_Floquet : 
The code calculates and saves Van der Pol's limit cycle, frequency, floquet vectors, etc. 
The calculation is in section1, drawing figures and saving the figures and data is in section2. The figure is corresponding to fig2 in the article*.

VAN_simple : It is necessary to run VAN_Floquet in advance.
The code simulates Zlotnik optimal inputs(2013) for van der Pol model. 
The first section is calculating lagrange multipliers for inputs. 
The second section is simulating dynamics with the input and holding simulation results. 
The third section is calculating phase coupling function and saving simulation data. 

VAN_feedback : It is necessary to run VAN_Floquet in advance.
The code simulates amplitude-feedback method for van der Pol model. 
The structure of section is the same as "VAN_simple"

VAN_penalty : It is necessary to run VAN_Floquet in advance.
The code simulates amplitude-penalty method for van der Pol model. 
The structure of section is the same as "VAN_simple"

VAN_tangent : It is necessary to run VAN_Floquet in advance.
The code simulates amplitude-penalty method for van der Pol model. 
The structure of section is the same as "VAN_simple"

VAN_plot : It is necessary to save simulation data for drawing figures in advance.
The code plots simulation results for van der Pol model. 
The first section is setting matplotlib parameters. 
The second section is drawing and saving simulation results of the amplitude-feedback method. The figure is corresponding to fig3 in the article*.
The third section is drawing and saving simulation results of the amplitude-penalty method. The figure is corresponding to fig3 in the article*.
The fourth section is drawing and saving simulation results of the amplitude-feedback and amplitude-penalty method. The figure is corresponding to fig4 in the article*.
The sixth section is drawing and saving simulation results of the tangential periodic input. The figure is corresponding to fig8 in the article*.

calculate_arnold_tongue : It is necessary to run VAN_Floquet in advance.
The code simulates three method, drawing arnold tongues and saving the figures and data. 
The first section is setting three functions for next sections. 
The second section is simulating Zlotnik optimal inputs(2013), drawing arnold tongues and saving the figures and data. The figure is corresponding to fig5 in the article*.
The third section is simulating amplitude-feedback method and so.
The fourth section is simulating amplitude-penalty method and so.
The fifth section is drawing and saving colorbar.


# Willamowski-Rossler

The main code is "WR_Floquet", "WR_simple", "WR_feedback", "WR_penalty", and "WR_plot". 
The structure of these code are the same as the above.
The figure is corresponding to fig6 and fig7 in the article*.  


# Fitz Hugh nagumo

The main code is "FHN_Floquet", "FHN_simple", "FHN_feedback", and "FHN_penalty". 
The structure of these code are the same as the above.
This model is not attached in the article. 


*Takata, Shohei, Yuzuru Kato, and Hiroya Nakao. "Fast optimal entrainment of limit-cycle oscillators by strong periodic inputs via phase-amplitude reduction and Floquet theory." arXiv preprint arXiv:2104.09944 (2021).


###############
## Questions ##
############### 

If you have any questions, don't hesitate to email takata.s.ae@m.titech.ac.jp
