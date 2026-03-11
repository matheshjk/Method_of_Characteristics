# Method of Characteristics Solver for Underexpanded Jet

## Overview

This project implements a **Python-based Method of Characteristics (MoC)
solver** to compute the two-dimensional flow field of an **underexpanded
supersonic jet**.

The solver models the expansion of a jet exiting a nozzle where the exit
pressure is higher than ambient pressure. Expansion waves form at the
nozzle lip and propagate through the flow, reflecting across the jet
axis until shock formation occurs.

The code computes:

-   Characteristic line intersections
-   Flow deflection angles
-   Mach number distribution
-   Pressure variation
-   Mach contour visualization

The resulting solution provides a physically consistent representation
of the **expansion fan and reflection structure** in the jet plume.

------------------------------------------------------------------------

# Physics Background

When a supersonic jet exits a nozzle with **exit pressure greater than
ambient pressure**, the flow undergoes **isentropic expansion** through
a Prandtl--Meyer expansion fan.

This produces two families of characteristic lines:

-   **Γ⁺ (positive characteristics)**
-   **Γ⁻ (negative characteristics)**

Along these characteristics, the compatibility equations apply:

ν + φ = constant (Γ⁻)

ν − φ = constant (Γ⁺)

Where:

-   ν = Prandtl--Meyer function
-   φ = flow turning angle

The solution becomes invalid when characteristics intersect, indicating
**shock formation in the jet plume**.

------------------------------------------------------------------------

# Numerical Method

The solver performs the following steps.

## 1. Initial Expansion Fan

An expansion fan is generated from the nozzle exit.

The code discretizes the fan using a user-defined number of
characteristics.

For each characteristic:

-   Prandtl--Meyer angles are computed
-   Mach numbers are obtained via **bisection inversion of the PM
    function**
-   Characteristic slopes are calculated

------------------------------------------------------------------------

## 2. Characteristic Intersections

The intersections between Γ⁺ and Γ⁻ lines are computed by solving the
equations of two straight lines.

The script uses **SymPy** to compute the intersection coordinates.

------------------------------------------------------------------------

## 3. Reflection Regions

The jet expansion includes successive reflection regions:

1.  Initial Γ⁻ expansion fan
2.  Γ⁺ reflection region
3.  Secondary Γ⁻ reflections

Each region is solved using compatibility equations along the
characteristics.

------------------------------------------------------------------------

## 4. Flow Property Computation

For every node in the characteristic mesh, the solver computes:

-   Flow turning angle
-   Prandtl--Meyer angle
-   Mach number
-   Pressure via isentropic relations
-   Mach angle

------------------------------------------------------------------------

## 5. Mach Contour Visualization

The flow field is visualized using triangular interpolation.

The solver generates:

-   Characteristic line network
-   Mach number contour plot
-   Jet boundary

------------------------------------------------------------------------

# Example Simulation Parameters

Default parameters in the script:

  Parameter          Value
  ------------------ ---------------
  Exit Mach number   2
  Exit pressure      1.2 × ambient
  Gas                Air
  γ                  1.4

The number of characteristics is chosen by the user at runtime.

------------------------------------------------------------------------

# Running the Code

Install required libraries:

    pip install numpy matplotlib sympy

Run the script:

    python Task_2_MachPlot.py

The program will prompt:

    Number of characteristics:

Higher numbers increase resolution of the characteristic mesh.

------------------------------------------------------------------------

# Example Output

The solver produces a **Mach number contour map** of the jet plume
showing:

-   Prandtl--Meyer expansion fan
-   Characteristic reflections
-   Increasing Mach number along the jet axis
-   Regions where the MoC solution breaks down due to shock formation

------------------------------------------------------------------------

# Key Skills Demonstrated

-   Compressible flow modelling
-   Method of Characteristics
-   Gas dynamics
-   Numerical root finding
-   Python scientific computing
-   Flow visualization with Matplotlib

------------------------------------------------------------------------

# Author

Mathesh Babu JK\
MSc Aerospace Engineering -- TU Delft\
Aerodynamics \| CFD \| Compressible Flow
****
