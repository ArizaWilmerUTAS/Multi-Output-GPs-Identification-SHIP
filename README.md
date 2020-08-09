# System indetification of Ship Dynamic System with Multi Output Gaussian Processes

This code is based in the work of Prof. Mauricio A. Alvarez and Prof. Neil D. Lawrence [SheffieldML](https://github.com/SheffieldML/multigp) (some modification were made such that Gaussian Processes with dimension greater than 2 can be run) and the [work](https://github.com/Dynamic-Systems-and-GP/GPdyn) of Prof. Dr. Juš Kocijan. This is an implementation of dynamic system identification with multi-output Gaussian Processes. Both dependencies need to be installed before running the code. An Underactuated Ship model was used to generate data for system identification. The system is a two input , four output system. A NARX architecture was used for the Multi-output Gaussian Processes as for the comparative Neural Networks.

If you want to cite this work please cite:

Ariza R., W., Leong, Z.Q., Nguyen, H. and Jayasinghe, S.G., 2018. Non-parametric dynamic system identification of ships using multi-output Gaussian Processes. Ocean Engineering, 166, pp.26-36.

<center><img src="figure9.png" width ="75%"><img src="figure10.png" width ="75%"><img src="figure11.png" width ="75%"><img src="figure12.png" width ="75%">
<br>Prediction from Multi_output GPs by algorithm of Naive Simulation with full data compared to mathematical model, a) controlled surge acceleration, b) induced sway speed, c) controlled yaw speed, and d) induced roll speed </center>


## Getting Started

To run the System idetification, please run the file: MultiOutputGPSSI.m





## Author


Wilmer Ariza Ramirez

Australian Maritime College, 
University of Tasmania, Newnham TAS 7248, Australia

Wilmer.ArizaRamirez@utas.edu.au 

 

