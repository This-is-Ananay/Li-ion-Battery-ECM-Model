# Li-ion-Battery-ECM-Model
This project focuses on the development and implementation of a lithium-ion battery 1-RC Thevenin equivalent circuit model using MATLAB. The main objective is to characterize behaviour, estimate internal parameters and validate the voltage response on multiple dynamic drive cycle tests. It features Pseudo OCV-SOC relationship, parameters estimation, static hystersis effect and validation.

## MODEL DEVELOPMENT
The modelling workflow consists of the following stages:

1- Experimental battery data processing.
2- Pseudo OCV-SOC characterization.
3- Battery parameter identification.
4- Equivalent circuit model implementation.
5- Dynamic drive cycle validation.

## Repository Structure

- [`data/`](Data/) - Experimental battery datasets
- [`model/`](model/) - Battery model implementation
- [`test/`](test/) - Model validation scripts
- [`lookup/`](Lookup/) - Lookup tables
- [`results/`](plots/) - Simulation plots
- [`run_all.m`](Run_all.m) - Main execution script

## DATASET
This dataset presents the experimental campaign for a batch of eleven prismatic CALB L148N58A lithium-ion B-grade battery cells with a nominal capacity of 58 Ah. Particulary cell 59294 is used for the project with a testing temperature of 25°C.
Charging and discharging test are each done at C/20 and HPPC at 1C. 

Link: https://data.mendeley.com/datasets/ycx459r5c3/2

## IMPLEMENTATION

## PSEUDO OCV-SOC CHARACTERSTICS
The dataset is not a perfect dataset for OCV-SOC characteristics as the no proper rest is given during testing for terminal voltage to fall near OCV. That's why a pseudo ocv-soc characterstics are obtained assuming at low rate C/20 terminal voltage is near equivalent to ocv.

Firstly two separate look-up are created for Charging and discharging with respect to terminal voltages. Then both are interpolated and the average of terminal voltages is taken to formulate OCV-SOC characteristics. A static hysteresis term(not corrected) is also added which at this stage is equivalent to half gap between charge and discharge curves. This will be corrected later on.

<p align="center">
  <img src="Screenshots/simulation/Schematics.png" width="48%">
  <img src="Screenshots/simulation/WITHOUT_CLAMP.png" width="48%">
</p>





## RESULTS





