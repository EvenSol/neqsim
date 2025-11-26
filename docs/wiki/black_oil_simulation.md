# Black-oil simulations in NeqSim

This page outlines how to build and run black-oil style simulations in NeqSim. It covers how to supply PVT data, create streams, execute flash calculations, and connect the models to process equipment.

## Core building blocks

NeqSim models black-oil systems with a few focused classes:

- **`BlackOilPVTTable`** stores pressure-indexed records with Rs, Rv, Bo/Bg/Bw and viscosities, and linearly interpolates values across pressures.【F:src/main/java/neqsim/blackoil/BlackOilPVTTable.java†L8-L210】
- **`BlackOilFlash`** performs the black-oil flash at a given pressure/temperature using the PVT table and reference standard-condition densities.【F:src/main/java/neqsim/blackoil/BlackOilFlash.java†L3-L95】
- **`SystemBlackOil`** is a lightweight stream representation that tracks standard totals, pressure, and temperature, and exposes convenience getters for PVT and flash-derived properties.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L3-L313】
- **`BlackOilConverter`** builds a `BlackOilPVTTable` and matching `SystemBlackOil` from a compositional (EOS-based) `SystemInterface` across a pressure grid.【F:src/main/java/neqsim/blackoil/BlackOilConverter.java†L11-L196】
- **`EclipseBlackOilImporter`** parses ECLIPSE decks (PVTO/PVTG/PVTW blocks plus densities) to create a PVT table and initialized black-oil system, handling METRIC/FIELD/LAB units.【F:src/main/java/neqsim/blackoil/io/EclipseBlackOilImporter.java†L17-L160】【F:src/main/java/neqsim/blackoil/io/EclipseBlackOilImporter.java†L340-L441】
- **`BlackOilSeparator`** provides a simple equilibrium separator that splits an inlet `SystemBlackOil` into oil, gas, and water outlets at a target pressure and temperature.【F:src/main/java/neqsim/process/equipment/blackoil/BlackOilSeparator.java†L3-L130】

## Supplying PVT data

### Importing an ECLIPSE deck
Use the importer when you have PVTO/PVTG/PVTW blocks and standard-condition densities:

```java
Path deck = Paths.get("/path/to/blackoil.INC");
EclipseBlackOilImporter.Result res = EclipseBlackOilImporter.fromFile(deck);
SystemBlackOil blackOil = res.system; // already seeded with bubble-point pressure and Rs
```
The importer converts FIELD units, interpolates gas/water properties where needed, and logs what it found in `res.log`. The resulting `res.pvt` holds the interpolated table and `res.bubblePoint` identifies the bubble-point pressure.【F:src/main/java/neqsim/blackoil/io/EclipseBlackOilImporter.java†L17-L160】【F:src/main/java/neqsim/blackoil/io/EclipseBlackOilImporter.java†L340-L441】

### Converting a compositional fluid
If you start from a standard NeqSim `SystemInterface`, you can generate a black-oil table over a pressure grid:

```java
double[] pGrid = new double[] {50, 100, 200, 300};
BlackOilConverter.Result res = BlackOilConverter.convert(eosFluid, 350.0, pGrid, 1.01325, 288.15);
SystemBlackOil blackOil = res.blackOilSystem; // totals set from stock-tank flash
```
The converter performs TP flashes at each pressure to compute Rs/Rv, formation volume factors, viscosities, and estimates standard-condition densities when possible. Values above the bubble point reuse saturated Rs and the last gas properties to remain consistent with black-oil conventions.【F:src/main/java/neqsim/blackoil/BlackOilConverter.java†L11-L196】

### Creating a table manually
You can also build a table directly if you already have PVT data points:

```java
List<BlackOilPVTTable.Record> recs = List.of(
    new BlackOilPVTTable.Record(100.0, 90.0, 1.25, 1.2e-3, 0.005, 1.0e-5, 0.0, 1.02, 0.5e-3),
    new BlackOilPVTTable.Record(200.0, 80.0, 1.20, 1.3e-3, 0.004, 1.0e-5, 0.0, 1.01, 0.5e-3)
);
BlackOilPVTTable pvt = new BlackOilPVTTable(recs, 150.0);
SystemBlackOil stream = new SystemBlackOil(pvt, 820.0, 1.1, 995.0);
```
`BlackOilPVTTable` linearly interpolates between pressures and clamps below/above the supplied range. When pressure exceeds the bubble point, `RsEffective` keeps Rs fixed at the bubble-point value as expected for black-oil calculations.【F:src/main/java/neqsim/blackoil/BlackOilPVTTable.java†L88-L210】

## Running flash calculations

1. Set the thermodynamic state on the stream:
   ```java
   stream.setPressure(120.0);
   stream.setTemperature(360.0);
   stream.setStdTotals(1.0, 100.0, 0.05); // Sm3 oil, Sm3 gas, Sm3 water
   ```
2. Call `stream.flash()` to compute phase split and property results. The flash returns reservoir volumes (V_o/V_g/V_w), viscosities, densities, and current Rs/Rv/Bo/Bg/Bw.
3. Use the convenience getters for quick access (e.g., `getBo()`, `getOilDensity()`, `getGasViscosity()`, `getOilReservoirVolume()`). Subsequent getter calls reuse the cached flash result until you change pressure, temperature, or totals.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L52-L313】

`BlackOilFlash` applies the standard black-oil material balance: it determines free gas and remaining solution-oil volumes from the standard totals, scales them with Bo/Bg/Bw, and computes phase densities by combining standard densities with Rs/Rv.【F:src/main/java/neqsim/blackoil/BlackOilFlash.java†L32-L95】

## Using black-oil equipment

For quick flow-sheet style simulations, attach a `SystemBlackOil` to `BlackOilSeparator`:

```java
BlackOilSeparator sep = new BlackOilSeparator("HP separator", stream, 70.0, 330.0);
sep.run();
SystemBlackOil oil = sep.getOilOut();
SystemBlackOil gas = sep.getGasOut();
SystemBlackOil water = sep.getWaterOut();
```
`run()` performs a black-oil flash at the outlet conditions, then distributes standard volumes to oil, gas, and water outlet streams, honoring Rs and Rv for dissolved phases.【F:src/main/java/neqsim/process/equipment/blackoil/BlackOilSeparator.java†L101-L130】

## Method reference

- **`SystemBlackOil`**
  - `setStdTotals(double oilSm3, double gasSm3, double waterSm3)` – store standard-condition totals and invalidate cached flash results.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L76-L90】
  - `setPressure(double P)` / `setTemperature(double T)` – update state and clear cached flash results.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L52-L74】
  - `flash()` – execute the black-oil flash using current state and totals; repeated calls reuse the last result when state is unchanged.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L147-L159】
  - `getBo()`, `getBg()`, `getBw()`, `getRs()`, `getRv()` – interpolate PVT properties at the current pressure.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L161-L214】
  - `getOilDensity()`, `getGasDensity()`, `getWaterDensity()` and `get*Viscosity()` – access flash-computed densities and viscosities for each phase.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L220-L280】
  - `getOilReservoirVolume()`, `getGasReservoirVolume()`, `getWaterReservoirVolume()` – reservoir volumes after the flash.【F:src/main/java/neqsim/blackoil/SystemBlackOil.java†L282-L313】
- **`BlackOilPVTTable`**
  - `Record` constructor captures a PVT point; interpolation methods `Rs(p)`, `Bo(p)`, `Bg(p)`, `Bw(p)`, `Rv(p)`, and `mu_*` apply linear interpolation; `RsEffective(p)` holds Rs constant above the bubble point.【F:src/main/java/neqsim/blackoil/BlackOilPVTTable.java†L15-L210】
- **`BlackOilFlash`**
  - `flash(P, T, Otot_std, Gtot_std, W_std)` – computes free gas, remaining oil, reservoir volumes, viscosities, densities, and returns a `BlackOilFlashResult` DTO populated with Bo/Bg/Bw and Rs/Rv at the specified pressure.【F:src/main/java/neqsim/blackoil/BlackOilFlash.java†L32-L95】
- **Integration helpers**
  - `EclipseBlackOilImporter.fromFile(Path)` or `fromReader(Reader)` parse deck data, convert units when needed, assemble a `BlackOilPVTTable`, and seed a `SystemBlackOil` with bubble-point conditions and default totals.【F:src/main/java/neqsim/blackoil/io/EclipseBlackOilImporter.java†L17-L160】【F:src/main/java/neqsim/blackoil/io/EclipseBlackOilImporter.java†L340-L441】
  - `BlackOilConverter.convert(...)` derives a black-oil representation from any compositional NeqSim fluid by flashing across a pressure grid and normalizing Rs above the bubble point.【F:src/main/java/neqsim/blackoil/BlackOilConverter.java†L30-L136】
