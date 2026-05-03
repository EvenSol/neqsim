---
title: JSON Format for ProcessSystem and ProcessModel
description: Accurate reference for NeqSim JSON process builder input, validation, and large flowsheet setup with mixers, splitters, and recycles.
---

# JSON Format for ProcessSystem and ProcessModel

This page documents the **actual JSON format** consumed by `ProcessSystem.fromJson(...)` and
`ProcessSystem.fromJsonAndRun(...)`.

## Quick answer: is there a JSON verifier?

Yes:

- `ProcessSystem.validateJson(String json)`
- `ProcessJsonValidator.validate(String json)`

Use validation as a pre-flight check before calling `fromJson()` or `fromJsonAndRun()`.

## 1. Root JSON structure (builder input)

The builder expects a root object with these common keys:

- `fluid` (optional): single default fluid definition
- `fluids` (optional): named fluids map (for `fluidRef` on streams)
- `process` (**required**): array of unit definitions
- `autoRun` (optional): `true` to run immediately

Minimal valid structure:

```json
{
  "fluid": {
    "model": "SRK",
    "temperature": 298.15,
    "pressure": 50.0,
    "mixingRule": "classic",
    "components": {
      "methane": 0.85,
      "ethane": 0.10,
      "propane": 0.05
    }
  },
  "process": [
    {
      "type": "Stream",
      "name": "feed",
      "properties": {
        "flowRate": [50000.0, "kg/hr"]
      }
    }
  ]
}
```

> Important: current builder input uses `process` + `properties` + `inlet`/`inlets` patterns.

## 2. Unit definition format

Each object in `process` should include:

- `type` (required)
- `name` (strongly recommended and should be unique)
- wiring fields depending on unit type:
  - `inlet`: single reference string
  - `inlets`: array of reference strings
- `properties` object for operating conditions and setpoints

Example:

```json
{
  "type": "Compressor",
  "name": "Comp",
  "inlet": "HP Sep.gasOut",
  "properties": {
    "outletPressure": [80.0, "bara"]
  }
}
```

## 3. Stream addressing and ports

References are by unit name, optionally with port suffix:

- `"feed"` (stream unit)
- `"HP Sep.gasOut"`
- `"HP Sep.liquidOut"`
- `"Comp.out"` / `"Comp.outStream"`

The resolver trims surrounding whitespace.

## 4. Large process JSON example (splitter + mixer + recycle)

```json
{
  "fluid": {
    "model": "SRK",
    "temperature": 298.15,
    "pressure": 70.0,
    "mixingRule": "classic",
    "components": {
      "methane": 0.88,
      "ethane": 0.08,
      "propane": 0.04
    }
  },
  "autoRun": true,
  "process": [
    {"type": "Stream", "name": "feed", "properties": {"flowRate": [250000.0, "kg/hr"]}},

    {"type": "Splitter", "name": "Feed Splitter", "inlet": "feed",
      "properties": {"splitFactors": [0.60, 0.40]}},

    {"type": "Cooler", "name": "Train A Cooler", "inlet": "Feed Splitter.splitStream_0",
      "properties": {"outTemperature": [20.0, "C"]}},

    {"type": "ThrottlingValve", "name": "Train B Valve", "inlet": "Feed Splitter.splitStream_1",
      "properties": {"outletPressure": [65.0, "bara"]}},

    {"type": "Mixer", "name": "Combined Mixer",
      "inlets": ["Train A Cooler.out", "Train B Valve.out"]},

    {"type": "Separator", "name": "HP Sep", "inlet": "Combined Mixer.out"},

    {"type": "Compressor", "name": "Recycle Comp", "inlet": "HP Sep.gasOut",
      "properties": {"outletPressure": [72.0, "bara"]}},

    {"type": "Recycle", "name": "Gas Recycle", "inlet": "Recycle Comp.out"},

    {"type": "Mixer", "name": "Feed + Recycle",
      "inlets": ["feed", "Gas Recycle.out"]}
  ]
}
```

## 5. Validation workflow (recommended)

Java:

```java
String json = ...;
ProcessJsonValidator.ValidationReport report = ProcessSystem.validateJson(json);
if (!report.isValid()) {
  throw new IllegalArgumentException("Invalid process JSON: " + report.getErrors());
}
SimulationResult result = ProcessSystem.fromJsonAndRun(json);
```

Python (via JPype wrapper):

```python
report = ns.ProcessSystem.validateJson(json_text)
if not report.isValid():
    raise ValueError(str(report.getErrors()))
result = ns.ProcessSystem.fromJsonAndRun(json_text)
```

## 6. What validator checks today

Errors:
- invalid JSON
- missing `process` array
- missing `type`
- missing `name`
- duplicate names

Warnings:
- inlet/inlets/streams input references that do not resolve to known unit/port names

## 7. ProcessModel lifecycle JSON

For lifecycle snapshots/versioning of full `ProcessModel` and `ProcessSystem` states,
see:

- `docs/process/lifecycle/process_model_lifecycle.md`
- `docs/simulation/process_serialization.md`

Builder input JSON (`fromJson`) and lifecycle snapshot JSON (`ProcessModelState`) have different purposes.
