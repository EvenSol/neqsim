---
title: JSON Format for Process Systems and Process Models
description: Detailed schema-style guide for building NeqSim ProcessSystem and ProcessModel objects from JSON, including large flowsheets with recycles, mixers, and splitters.
---

# JSON Format for Process Systems and Process Models

This guide documents the **JSON** format used to build and execute NeqSim process simulations with:

- `ProcessSystem.fromJsonAndRun(...)` for a single flowsheet
- `ProcessModel` JSON state patterns for multi-area plants

It focuses on practical usage for large configurations, including **mixers**, **splitters**, and **recycle loops**.

## 1) Core Builder Pattern (`ProcessSystem.fromJsonAndRun`)

At minimum, the JSON has:

- `schemaVersion` (recommended)
- `units`: a list of equipment definitions
- optional `connections`: explicit wiring metadata

Typical Python execution:

```python
import json

ProcessSystem = ns.ProcessSystem
result = ProcessSystem.fromJsonAndRun(json.dumps(process_json))
if result.isSuccess():
    process = result.getProcessSystem()
```

Typical Java execution:

```java
SimulationResult result = ProcessSystem.fromJsonAndRun(jsonString);
if (result.isSuccess()) {
    ProcessSystem process = result.getProcessSystem();
}
```

## 2) Unit Object Structure

Each unit in `units` generally follows:

```json
{
  "name": "HP Separator",
  "type": "Separator",
  "inputs": {
    "pressure": {"value": 85.0, "unit": "bara"}
  },
  "streams": {
    "inlet": "Feed",
    "gasOut": "HP Separator.gasOutStream",
    "liquidOut": "HP Separator.liquidOutStream"
  }
}
```

Conventions:

- `name`: unique tag used in references.
- `type`: NeqSim equipment type/factory key.
- `inputs`: input variables (setpoints, targets, specs).
- `streams`: named stream references for in/out wiring.

## 3) Stream References and Addressing

NeqSim supports dot-addressing for named outlet streams:

- `"HP Separator.gasOutStream"`
- `"Compressor.outStream"`

For robust automation and debugging:

- keep unit names stable and human-readable
- avoid ambiguous aliases
- prefer explicit stream names in larger flowsheets

## 4) Large Process Example (Mixers, Splitters, Recycle)

The example below shows a realistic topology:

- feed split into two trains (`Splitter`)
- train recombination (`Mixer`)
- recycle gas returned upstream (`Recycle`)
- final export compression

```json
{
  "schemaVersion": "1.0",
  "name": "Large Gas Process with Recycle",
  "units": [
    {
      "name": "Feed",
      "type": "Stream",
      "fluid": {
        "temperature": {"value": 35.0, "unit": "C"},
        "pressure": {"value": 90.0, "unit": "bara"},
        "components": [
          {"name": "methane", "amount": 0.88},
          {"name": "ethane", "amount": 0.08},
          {"name": "propane", "amount": 0.04}
        ],
        "mixingRule": "classic"
      },
      "inputs": {
        "flowRate": {"value": 250000.0, "unit": "kg/hr"}
      }
    },
    {
      "name": "Feed Splitter",
      "type": "Splitter",
      "inputs": {
        "splitFactors": [0.60, 0.40]
      },
      "streams": {
        "inlet": "Feed",
        "splitOut1": "Feed Splitter.splitStream_0",
        "splitOut2": "Feed Splitter.splitStream_1"
      }
    },
    {
      "name": "Train A Cooler",
      "type": "Cooler",
      "inputs": {
        "outTemperature": {"value": 20.0, "unit": "C"}
      },
      "streams": {
        "inlet": "Feed Splitter.splitStream_0"
      }
    },
    {
      "name": "Train B Valve",
      "type": "ThrottlingValve",
      "inputs": {
        "outletPressure": {"value": 65.0, "unit": "bara"}
      },
      "streams": {
        "inlet": "Feed Splitter.splitStream_1"
      }
    },
    {
      "name": "Combined Mixer",
      "type": "Mixer",
      "streams": {
        "inlets": [
          "Train A Cooler.outStream",
          "Train B Valve.outStream"
        ],
        "outlet": "Combined Mixer.outStream"
      }
    },
    {
      "name": "HP Separator",
      "type": "Separator",
      "streams": {
        "inlet": "Combined Mixer.outStream",
        "gasOut": "HP Separator.gasOutStream",
        "liquidOut": "HP Separator.liquidOutStream"
      }
    },
    {
      "name": "Recycle Gas Compressor",
      "type": "Compressor",
      "inputs": {
        "outletPressure": {"value": 92.0, "unit": "bara"}
      },
      "streams": {
        "inlet": "HP Separator.gasOutStream"
      }
    },
    {
      "name": "Gas Recycle",
      "type": "Recycle",
      "streams": {
        "inlet": "Recycle Gas Compressor.outStream",
        "outlet": "Recycle Seed"
      }
    },
    {
      "name": "Recycle Seed Mixer",
      "type": "Mixer",
      "streams": {
        "inlets": ["Feed", "Recycle Seed"],
        "outlet": "Feed With Recycle"
      }
    },
    {
      "name": "Export Compressor",
      "type": "Compressor",
      "inputs": {
        "outletPressure": {"value": 120.0, "unit": "bara"}
      },
      "streams": {
        "inlet": "HP Separator.gasOutStream"
      }
    }
  ],
  "connections": [
    {"from": "Feed", "to": "Feed Splitter", "type": "MATERIAL", "label": "Main Feed"},
    {"from": "Recycle Seed Mixer", "to": "Feed Splitter", "type": "MATERIAL", "label": "Feed+Recycle"},
    {"from": "Feed Splitter", "to": "Train A Cooler", "type": "MATERIAL", "label": "Train A"},
    {"from": "Feed Splitter", "to": "Train B Valve", "type": "MATERIAL", "label": "Train B"},
    {"from": "Train A Cooler", "to": "Combined Mixer", "type": "MATERIAL", "label": "A to Mixer"},
    {"from": "Train B Valve", "to": "Combined Mixer", "type": "MATERIAL", "label": "B to Mixer"},
    {"from": "Combined Mixer", "to": "HP Separator", "type": "MATERIAL", "label": "To Separation"},
    {"from": "HP Separator", "to": "Recycle Gas Compressor", "type": "MATERIAL", "label": "Gas Recycle Path"},
    {"from": "Recycle Gas Compressor", "to": "Gas Recycle", "type": "MATERIAL", "label": "Recycle Closure"},
    {"from": "Gas Recycle", "to": "Recycle Seed Mixer", "type": "MATERIAL", "label": "Back to Front"}
  ]
}
```

## 5) Recycle Guidance

For recycle convergence:

- Provide a physically reasonable recycle seed stream early in model setup.
- Keep pressure levels consistent around the loop.
- Add one recycle loop at a time when debugging large models.

## 6) ProcessModel JSON State (Multi-Area Plants)

For multiple areas (e.g., separation + compression), use `ProcessModelState` JSON snapshots:

- top-level metadata (`name`, `version`, `createdAt`)
- per-area process states (`processStates`)
- optional inter-area links (`interProcessConnections`)
- execution settings (`executionConfig`)

This format is intended for lifecycle/versioning and can be compared with `ProcessModelState.compare(...)`.

## 7) Validation Checklist

Before committing a large JSON model:

1. All unit `name` values are unique.
2. Every referenced stream address exists.
3. Fluid setup includes `mixingRule` (`"classic"` recommended baseline).
4. Recycle loop has a valid seed and pressure-feasible return path.
5. Simulation runs cleanly with no unresolved connections.

## 8) Related Documentation

- `docs/process/lifecycle/process_model_lifecycle.md`
- `docs/process/README.md`
- `docs/process/controllers.md`
- `docs/process/equipment/mixers_splitters.md`
- `docs/process/equipment/util/recycles.md`
