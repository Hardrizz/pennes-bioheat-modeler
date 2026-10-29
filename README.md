# pennes-bioheat-modeler
A biomedical engineering numerical modeling tool that simulates 1D heat transfer in living biological tissues using the Pennes Bioheat Equation to analyze temperature distributions during thermal therapies.
Tissue Thermal Profile & Pennes Bioheat Modeler

## Description
A biomedical engineering numerical solver designed to model heat transfer within living biological tissues. This application implements the 1D steady-state *Pennes Bioheat Equation*, accounting for thermal conduction, tissue metabolism, and capillary blood perfusion to simulate temperature profiles during hyperthermia treatments or thermal ablation therapies.

## Key Features
* *Finite Difference Solver:* Numerically solves differential equations governing heat propagation across biological barriers.
* *Perfusion & Metabolism Modeling:* Incorporates terms for blood perfusion cooling/warming ($w_b c_b [T_a - T]$) and metabolic heat generation ($Q_{met}$).
* *Thermal Threshold Tracking:* Evaluates temperature gradients relative to critical core body standards and clinical ablation benchmarks.
* *Visual Profile Generation:* Exports high-resolution graphs illustrating spatial temperature distribution across tissue depth.

## Mathematical Formulation
The steady-state 1D Pennes Bioheat Equation solved:
$$k \frac{d^2 T}{dx^2} - w_b c_b (T - T_a) + Q_{met} = 0$$

Where:
* $k$ = Thermal conductivity of tissue ($W/m \cdot K$)
* $T$ = Tissue temperature ($^\circ C$)
* $w_b$ = Blood perfusion rate ($kg/m^3 \cdot s$)
* $c_b$ = Specific heat of blood ($J/kg \cdot K$)
* $T_a$ = Arterial blood temperature ($^\circ C$)
* $Q_{met}$ = Metabolic heat source ($W/m^3$)

