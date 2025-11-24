# TurboExpander–Compressor modelling in NeqSim

This note summarizes how the coupled turboexpander–compressor unit (`TurboExpanderCompressor`) is
modelled in NeqSim and how to use it from Python. Equations are expressed in the same symbols as the
implementation to ease cross-referencing.

## Energy balance and speed matching

The `run` method iterates on shaft speed `N` so that expander power equals compressor power plus
mechanical losses. For each Newton–Raphson iteration:

- The expander outlet state is obtained from an isentropic flash at the target outlet pressure
  (`PSflash`), giving the specific isentropic enthalpy drop
  \(h_s = (h_{\text{in}} - h_{\text{out}}^{s})\times 1000\;[\text{J/kg}]\).【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L182-L214】
- The blade tip speed \(U = \pi D N / 60\) and ideal spouting velocity \(C = \sqrt{2 h_s}\) define
  a non-dimensional velocity ratio \(u_c = U/(C\,u_{c,\text{design}})\). The expander isentropic
  efficiency factor from this velocity ratio is given by a parabola
  \(\eta_{u_c} = -3.56 (u_c-1)^2 + 1\).【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L192-L208】【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L468-L503】
- A volumetric-flow correction is optionally applied using the non-dimensional flow ratio
  \(q_n = \tfrac{Q\,60/N}{(Q/N)_{\text{design}}}\) and a spline-based efficiency multiplier
  \(\eta_{q_n}\) from the user-provided Q/N curve.【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L197-L224】【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L534-L586】
- The expander isentropic efficiency becomes \(\eta_s = \eta_{s,\text{design}}\,\eta_{u_c}\,\eta_{q_n}\)
  and expander power is \(W_{\text{exp}} = \dot m_{\text{exp}}\,h_s\,\eta_s\).【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L212-L224】
- Compressor polytropic head scales with speed and Q/N-dependent curve factor,
  \(H_p = H_{p,\text{design}}(N/N_{\text{design}})^2\,f_{H}(q_n)\), while polytropic efficiency is
  \(\eta_p = \eta_{p,\text{design}}\,f_{\eta}(q_n)\).【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L216-L224】【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L588-L669】
- Compressor power and bearing losses are computed as
  \(W_{\text{comp}} = \dot m_{\text{comp}} H_p /(\eta_p)\,1000\) and
  \(W_{\text{bear}} = P_{\text{bear,design}}(N/N_{\text{design}})^2\).【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L221-L225】
- Newton–Raphson updates the speed to drive the residual
  \(f(N) = W_{\text{exp}} - (W_{\text{comp}} + W_{\text{bear}})\) toward zero with a finite
  difference derivative. Iteration stops when the power mismatch is small or limits on iteration
  count and allowable speeds are reached.【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L226-L306】

After convergence the class back-calculates outlet stream states with an auxiliary `Expander` and
polytropic `Compressor` to make results available downstream in the flowsheet.【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L293-L306】

## Curve fits and flow-control logic

Efficiency and head multipliers can be provided as Q/N tables; the implementation sorts user points
and applies monotonic cubic Hermite interpolation with linear extrapolation outside the data range.
If no curve is supplied the factors default to 1.0, preserving design performance without scaling.
【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L505-L586】【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L588-L669】

The inlet-guide-vane (IGV) helper estimates required area from current mass flow, density, and half
of the stage enthalpy drop \(\Delta h_\text{IGV} = 0.5 h_s\). The velocity \(v = \sqrt{2\,\Delta h_\text{IGV}}\)
sets the area \(A = \dot m/(\rho v)\); the model caps opening between 0 and 1 and optionally enables
an enlarged maximum area using the configured expansion factor. The IGV state (opening, active area,
last stage drop) is cached for later queries.【F:src/main/java/neqsim/process/equipment/expander/TurboExpanderCompressor.java†L724-L773】

## Using the model from Python

The NeqSim Python package wraps the Java classes via JPype. A minimal script for running a
`TurboExpanderCompressor` looks like this:

```python
import jpype.imports
import neqsim
from neqsim.thermo.system import SystemSrkEos
from neqsim.thermodynamicoperations import ThermodynamicOperations
from neqsim.process.equipment.stream import Stream
from neqsim.process.equipment.expander import TurboExpanderCompressor
from neqsim.process.processmodel import ProcessSystem

# 1) Build the feed fluid and flash to the inlet condition
fluid = SystemSrkEos(298.15, 70.0)
fluid.addComponent("methane", 0.9)
fluid.addComponent("ethane", 0.1)
fluid.createDatabase()
ThermodynamicOperations(fluid).TPflash()

# 2) Assemble the process
process = ProcessSystem()
inlet = Stream("feed", fluid)
process.add(inlet)

tex = TurboExpanderCompressor("TEX", inlet)
tex.setExpanderOutPressure(40.0)         # bara
tex.setDesignSpeed(6850.0)               # rpm
tex.setCompressorDesignPolytropicHead(20.47)  # kJ/kg at design
tex.setCompressorDesignPolytropicEfficiency(0.81)
tex.setExpanderDesignIsentropicEfficiency(0.88)
process.add(tex)

# 3) Run and read results
process.run()
print("Matched speed (rpm):", tex.getSpeed())
print("Expander power (MW):", tex.getPowerExpander()/1e6)
print("Compressor discharge pressure (bara):",
      tex.getCompressorOutletStream().getPressure("bara"))
```

The JPype bootstrap performed by `import neqsim` starts the JVM with the packaged NeqSim JAR, so the
Java classes used above are available directly in Python. Replace composition, design parameters,
or performance curves to reproduce a specific machine map before calling `process.run()`.
