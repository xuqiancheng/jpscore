[![GitHub license](https://img.shields.io/badge/license-GPL-blue.svg)](https://raw.githubusercontent.com/JuPedSim/jpscore/master/LICENSE)

## Introduction

**jpscore** is the core module of [jupedsim](https://www.jupedsim.org/jupedsim_introduction.html) for preforming the simulation (i.e. computing the trajectories).  
This branch is used to demonstrate the work of paper [On the Effectiveness of the Measures in Supermarkets for Reducing Contact among Customers during COVID-19 Period](https://www.mdpi.com/2071-1050/12/22/9385).  
The implementation of this work is mainly in [VelocityModel.cpp](https://github.com/xuqiancheng/jpscore/blob/COVID19/libcore/src/math/VelocityModel.cpp).  
For more information about jupedsim, see [documentation](https://www.jupedsim.org/). 

## Installation

See [Build jupedsim from source](https://www.jupedsim.org/jupedsim_requirements.html) for reference.  
Contact q.xu@fz-juelich.de if you meet any problems.

## Quick start

See [Getting started with jupedsim](http://www.jupedsim.org/jpscore_introduction.html).

## Parameters setting

Please read this paper first: [On the Effectiveness of the Measures in Supermarkets for Reducing Contact among Customers during COVID-19 Period](https://www.mdpi.com/2071-1050/12/22/9385).  
Also the document if you want to run simulations in the confined room: [The model for simulations in confined rooms](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf).

For the parameters used in the simulation of this work: 

```
<covid_parameters covid_parameter_id="1">
	<infective>1</infective>
	<Model_param k="0.1" D="1.0" /> 
	<P mu="0.3" sigma="0.00" />
	<Q mu="0.3" sigma="0.00" />
	<alpha mu="1.0" sigma="0.00" />
	<stayTime mu="300" sigma="50"/>
</covid_parameters>
```
- infective: **0** means the agent is healthy, **1** means the agent is infected.

+ Model_param:
	- k: the parameter corresponding to ![](http://latex.codecogs.com/gif.latex?K) in [equation (3) of confined room case](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf). 
	- D: the parameter corresponding to ![](http://latex.codecogs.com/gif.latex?D) in [equation (2) of confined room case](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf).
- P:  the parameter corresponding to ![](http://latex.codecogs.com/gif.latex?P_j) in [equation (2) of confined room case](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf). 
- Q:  the parameter corresponding to ![](http://latex.codecogs.com/gif.latex?Q_i) in [equation (4) of confined room case](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf). 
- alpha:  the parameter corresponding to ![](http://latex.codecogs.com/gif.latex?\alpha_j) in [equation (2) of confined room case](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf) and  ![](http://latex.codecogs.com/gif.latex?\alpha_i) in [equation (3) of confined room case](https://github.com/xuqiancheng/jpscore/blob/COVID19/demos/scenario_1_covid19_confined_room/Confine_room_model.pdf). 
- stayTime:  the shopping time for customers in the supermarket scenario, which is the parameter corresponding to ![](http://latex.codecogs.com/gif.latex?t_i^\text{shop} ) in [equation (5) of supermarket case](https://www.mdpi.com/2071-1050/12/22/9385).  

```
<model_parameters>
	<covid>5</covid> 
	<contact_para number="50" checkTimeRate="0.1" checkDisRate="0.003" socialDis="0" />
</model_parameters>
```
- covid: **1** means run simulations in confined room, **2**  means run simulations in normal case, **5** means run simulations in supermarket case
+ contact_para:
	- number: the max allowable number of customers in the supermarket at the same time, corresponds to ![](http://latex.codecogs.com/gif.latex?L_\text{max}) in [supermarket case](https://www.mdpi.com/2071-1050/12/22/9385).  
	- checkTimeRate: the rate between checkout time and shopping time, ![](http://latex.codecogs.com/gif.latex?\alpha) in [equation (5) of supermarket case](https://www.mdpi.com/2071-1050/12/22/9385).
	- checkDisRate: the rate between disatance and the shopping time, ![](http://latex.codecogs.com/gif.latex?\beta) in [equation (7) of supermarket case](https://www.mdpi.com/2071-1050/12/22/9385). 
	- socialDis:  the social distance, ![](http://latex.codecogs.com/gif.latex?d_\text{rule}) in [equation (7) of supermarket case](https://www.mdpi.com/2071-1050/12/22/9385).

The meaning of other parameters can be found in [jpscore inifile](https://www.jupedsim.org/jpscore_inifile.html).  

## Demos

The following demos are provided for testing. They can be found in [demos](https://github.com/xuqiancheng/jpscore/tree/COVID19/demos).  


- scenario_1_covid19_confined_room  
  Infected and healthy agents move in a confined room, the probability of infection is calculated for each healthy agent.  
  <img src="demos\scenario_1_covid19_confined_room\confined_room.png" style="zoom: 25%;" />  
- scenario_2_covid19_bottleneck  
  Infected and healthy agents move through bottlenecks with different width under periodic boundary conditions, the probability of infection is calculated for each healthy agent.  
  <img src="demos\scenario_2_covid19_bottleneck\bottleneck.png"  />
+ scenario_3_covid19_supermarket  
  supermarket simulations for [On the Effectiveness of the Measures in Supermarkets for Reducing Contact among Customers during COVID-19 Period](https://www.mdpi.com/2071-1050/12/22/9385).  
  <img src="demos\scenario_3_covid19_supermarket\supermarket.png" style="zoom: 25%;" />  

## Support

We are heavily working on this project which means that:

- It's not done. We will be releasing new enhancements, bug fixes etc.
- We love your support. If you find any errors or have suggestions, please write an issue in our [issue-tracker](https://github.com/JuPedSim/jpscore/issues). We will try hard to fix it.
- Be patient. We are scientists and PhD/master students. Therefore, we primarily care about our research and theses.

If you meet any problems with this branch, feel free to contact q.xu@fz-juelich.de.

Enjoy!

