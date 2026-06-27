package neqsim.process.equipment.absorber;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.UUID;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import neqsim.process.equipment.stream.StreamInterface;
import neqsim.process.equipment.util.Recycle;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Dedicated amine solvent regenerator for CO2 solvent regeneration service.
 *
 * <p>
 * The model is a compact process-unit representation of a stripper/reboiler/condenser/water-knockout loop. It accepts a
 * rich amine feed and produces lean amine bottoms plus a wet overhead vapour. The calculation uses available NeqSim
 * flash operations to initialise stage thermodynamics and applies amine-regeneration design correlations for lean
 * loading, CO2 recovery, reboiler duty, condenser duty, steam demand, water loss and degradation warnings.
 * </p>
 *
 * <p>
 * The class is intended for flowsheet screening and absorber/regenerator recycle iteration. It exposes lean/rich loading
 * getters so {@link RateBasedAbsorber} and process modules can iterate to a consistent loading before a more rigorous
 * electrolyte column model is introduced.
 * </p>
 *
 * @author NeqSim
 */
public class AmineRegenerator extends SimpleAbsorber {
  private static final long serialVersionUID = 1000L;
  private static final Logger logger = LogManager.getLogger(AmineRegenerator.class);
  private static final double CO2_HEAT_OF_DESORPTION_KJ_PER_MOL = 72.0;
  private static final double WATER_LATENT_HEAT_KJ_PER_KG = 2257.0;
  private static final double STEAM_LATENT_HEAT_KJ_PER_KG = 2100.0;

  /** Regenerator specification mode. */
  public enum SpecificationMode {
    /** Solve to target lean loading. */
    LEAN_LOADING,
    /** Solve from available reboiler duty. */
    REBOILER_DUTY,
    /** Solve from target CO2 recovery. */
    CO2_RECOVERY,
    /** Solve from stripper pressure and stage count/height. */
    STRIPPER_PRESSURE,
    /** Solve from number of stages or packing height. */
    STAGES_OR_PACKING_HEIGHT
  }

  private StreamInterface richAmineFeed;
  private StreamInterface leanAmineBottoms;
  private StreamInterface overheadVapor;
  private StreamInterface condenserLiquid;
  private StreamInterface refluxStream;
  private StreamInterface waterKnockoutLiquid;
  private SpecificationMode specificationMode = SpecificationMode.LEAN_LOADING;
  private String amineType = "MDEA";
  private double richLoading = 0.45;
  private double targetLeanLoading = 0.05;
  private double leanLoading = 0.05;
  private double reboilerDuty = Double.NaN;
  private double calculatedReboilerDuty = 0.0;
  private double condenserDuty = 0.0;
  private double stripperPressure = 1.8;
  private double reboilerTemperatureC = 120.0;
  private double packingHeight = 12.0;
  private double co2Recovery = 0.90;
  private double overheadCO2Purity = 0.0;
  private double waterMakeupRate = 0.0;
  private double steamRate = 0.0;
  private double co2ReleasedMolesPerSec = 0.0;
  private double refluxRatio = 0.25;
  private boolean electrolyteFlashAttempted = false;
  private boolean electrolyteFlashSucceeded = false;
  private final List<String> degradationWarningFlags = new ArrayList<String>();

  /**
   * Creates an amine regenerator.
   *
   * @param name equipment name
   */
  public AmineRegenerator(String name) {
    super(name);
    setNumberOfStages(12);
  }

  /**
   * Creates an amine regenerator with rich amine feed.
   *
   * @param name equipment name
   * @param richAmineFeed rich amine feed stream
   */
  public AmineRegenerator(String name, StreamInterface richAmineFeed) {
    this(name);
    setRichAmineFeed(richAmineFeed);
  }

  /**
   * Sets rich amine feed and initializes outlet stream placeholders.
   *
   * @param stream rich amine stream
   */
  public void setRichAmineFeed(StreamInterface stream) {
    this.richAmineFeed = stream;
    this.leanAmineBottoms = stream.clone(getName() + " lean amine bottoms");
    this.overheadVapor = stream.clone(getName() + " overhead CO2/H2O vapor");
    this.condenserLiquid = stream.clone(getName() + " condenser liquid");
    this.refluxStream = stream.clone(getName() + " reflux water");
    this.waterKnockoutLiquid = stream.clone(getName() + " water knockout liquid");
  }

  public StreamInterface getRichAmineFeed() { return richAmineFeed; }
  public StreamInterface getLeanAmineBottoms() { return leanAmineBottoms; }
  public StreamInterface getOverheadVapor() { return overheadVapor; }
  public StreamInterface getCondenserLiquid() { return condenserLiquid; }
  public StreamInterface getRefluxStream() { return refluxStream; }
  public StreamInterface getWaterKnockoutLiquid() { return waterKnockoutLiquid; }
  public void setAmineType(String amineType) { this.amineType = amineType; }
  public String getAmineType() { return amineType; }
  public void setRichLoading(double richLoading) { this.richLoading = richLoading; }
  public double getRichLoading() { return richLoading; }
  public void setLeanLoading(double leanLoading) { this.targetLeanLoading = leanLoading; this.specificationMode = SpecificationMode.LEAN_LOADING; }
  public double getLeanLoading() { return leanLoading; }
  public void setReboilerDuty(double duty, String unit) { this.reboilerDuty = toKw(duty, unit); this.specificationMode = SpecificationMode.REBOILER_DUTY; }
  public double getReboilerDuty(String unit) { return fromKw(calculatedReboilerDuty, unit); }
  public double getCondenserDuty(String unit) { return fromKw(condenserDuty, unit); }
  public void setStripperPressure(double pressure, String unit) { this.stripperPressure = toBar(pressure, unit); this.specificationMode = SpecificationMode.STRIPPER_PRESSURE; }
  public double getStripperPressure(String unit) { return "Pa".equals(unit) ? stripperPressure * 1.0e5 : stripperPressure; }
  public void setPackingHeight(double packingHeight) { this.packingHeight = packingHeight; this.specificationMode = SpecificationMode.STAGES_OR_PACKING_HEIGHT; }
  public double getPackingHeight() { return packingHeight; }
  public void setCO2Recovery(double recovery) { this.co2Recovery = limit(recovery, 0.0, 0.999); this.specificationMode = SpecificationMode.CO2_RECOVERY; }
  public double getCO2Recovery() { return co2Recovery; }
  public double getOverheadCO2Purity() { return overheadCO2Purity; }
  public double getWaterMakeupRate(String unit) { return "kg/hr".equals(unit) ? waterMakeupRate * 3600.0 : waterMakeupRate; }
  public double getSteamRate(String unit) { return "kg/hr".equals(unit) ? steamRate * 3600.0 : steamRate; }
  public boolean isElectrolyteFlashAttempted() { return electrolyteFlashAttempted; }
  public boolean isElectrolyteFlashSucceeded() { return electrolyteFlashSucceeded; }
  public List<String> getDegradationWarningFlags() { return Collections.unmodifiableList(degradationWarningFlags); }

  /** {@inheritDoc} */
  @Override
  public void run(UUID id) {
    if (richAmineFeed == null) {
      logger.error("Rich amine feed must be set before running amine regenerator");
      return;
    }
    SystemInterface richSystem = richAmineFeed.getThermoSystem().clone();
    flashStage(richSystem);
    richLoading = estimateLoading(richSystem, richLoading);
    double amineMoles = Math.max(getComponentMoles(richSystem, amineType), 1.0e-12);
    co2Recovery = calculateRecoveryFromSpecification(richLoading, amineMoles);
    leanLoading = Math.max(0.0, richLoading * (1.0 - co2Recovery));
    co2ReleasedMolesPerSec = Math.max(0.0, (richLoading - leanLoading) * amineMoles);
    calculatedReboilerDuty = calculateReboilerDutyKw(co2ReleasedMolesPerSec);
    condenserDuty = calculateCondenserDutyKw(co2ReleasedMolesPerSec);
    waterMakeupRate = calculateWaterLossKgPerSec(co2ReleasedMolesPerSec);
    steamRate = calculatedReboilerDuty / STEAM_LATENT_HEAT_KJ_PER_KG;
    overheadCO2Purity = co2ReleasedMolesPerSec / Math.max(co2ReleasedMolesPerSec + waterMakeupRate / 0.018015, 1.0e-12);
    updateStreams(richSystem, id);
    updateWarnings();
    setCalculationIdentifier(id);
  }

  /**
   * Connects this regenerator lean outlet to a rate-based absorber and optionally receives absorber rich solvent.
   *
   * @param absorber absorber to connect
   */
  public void integrateWith(RateBasedAbsorber absorber) {
    absorber.setRegenerator(this);
    if (leanAmineBottoms != null) {
      absorber.addSolventInStream(leanAmineBottoms);
    }
  }

  /**
   * Adds a recycle from this regenerator lean outlet to an absorber solvent inlet stream.
   *
   * @param recycle recycle unit to configure
   */
  public void configureLeanAmineRecycle(Recycle recycle) {
    if (leanAmineBottoms != null) {
      recycle.addStream(leanAmineBottoms);
    }
  }

  private void flashStage(SystemInterface system) {
    electrolyteFlashAttempted = true;
    try {
      new ThermodynamicOperations(system).TPflash();
      system.initProperties();
      electrolyteFlashSucceeded = true;
    } catch (Exception ex) {
      electrolyteFlashSucceeded = false;
      logger.warn("Amine regenerator stage flash failed, continuing with correlation model: {}", ex.getMessage());
      system.init(1);
    }
  }

  private double calculateRecoveryFromSpecification(double inletLoading, double amineMoles) {
    if (specificationMode == SpecificationMode.REBOILER_DUTY && !Double.isNaN(reboilerDuty)) {
      double recoverable = reboilerDuty / Math.max(calculateReboilerDutyKw(inletLoading * amineMoles), 1.0e-12);
      return limit(recoverable, 0.0, 0.995);
    }
    if (specificationMode == SpecificationMode.CO2_RECOVERY) {
      return co2Recovery;
    }
    if (specificationMode == SpecificationMode.STRIPPER_PRESSURE || specificationMode == SpecificationMode.STAGES_OR_PACKING_HEIGHT) {
      double pressureFactor = limit(2.0 / Math.max(stripperPressure, 0.5), 0.60, 1.15);
      double stageFactor = 1.0 - Math.exp(-0.18 * Math.max(getNumberOfStages(), packingHeight / 1.0));
      return limit(0.985 * pressureFactor * stageFactor, 0.0, 0.995);
    }
    return limit(1.0 - targetLeanLoading / Math.max(inletLoading, 1.0e-12), 0.0, 0.995);
  }

  private double calculateReboilerDutyKw(double co2MolesPerSec) {
    double sensibleHeatKw = richAmineFeed.getFlowRate("kg/hr") / 3600.0 * 3.6 * Math.max(0.0, reboilerTemperatureC - 45.0) * 0.08;
    return co2MolesPerSec * CO2_HEAT_OF_DESORPTION_KJ_PER_MOL + waterMakeupRate * WATER_LATENT_HEAT_KJ_PER_KG
        + sensibleHeatKw;
  }

  private double calculateCondenserDutyKw(double co2MolesPerSec) {
    double waterKgPerSec = calculateWaterLossKgPerSec(co2MolesPerSec) * (1.0 + refluxRatio);
    return waterKgPerSec * WATER_LATENT_HEAT_KJ_PER_KG;
  }

  private double calculateWaterLossKgPerSec(double co2MolesPerSec) {
    return Math.max(0.0, co2MolesPerSec * 0.018015 * (1.0 + refluxRatio));
  }

  private void updateStreams(SystemInterface richSystem, UUID id) {
    SystemInterface leanSystem = richSystem.clone();
    removeComponentMoles(leanSystem, "CO2", Math.min(getComponentMoles(leanSystem, "CO2"), co2ReleasedMolesPerSec));
    leanSystem.setPressure(stripperPressure, "bara");
    leanSystem.setTemperature(273.15 + Math.min(reboilerTemperatureC, 135.0));
    leanSystem.init(1);
    leanAmineBottoms.setThermoSystem(leanSystem);
    leanAmineBottoms.setCalculationIdentifier(id);

    SystemInterface overheadSystem = richSystem.clone();
    for (int i = overheadSystem.getNumberOfComponents() - 1; i >= 0; i--) {
      String name = overheadSystem.getComponent(i).getName();
      double keep = "CO2".equals(name) ? co2ReleasedMolesPerSec : ("water".equals(name) ? waterMakeupRate / 0.018015 : 0.0);
      removeComponentMoles(overheadSystem, name, Math.max(0.0, getComponentMoles(overheadSystem, name) - keep));
    }
    overheadSystem.setPressure(stripperPressure, "bara");
    overheadSystem.setTemperature(273.15 + 100.0);
    overheadSystem.init(1);
    overheadVapor.setThermoSystem(overheadSystem);
    overheadVapor.setCalculationIdentifier(id);
    condenserLiquid.setThermoSystem(overheadSystem.clone());
    refluxStream.setThermoSystem(overheadSystem.clone());
    waterKnockoutLiquid.setThermoSystem(overheadSystem.clone());
  }

  private void updateWarnings() {
    degradationWarningFlags.clear();
    if (reboilerTemperatureC > maxRecommendedReboilerTemperature()) {
      degradationWarningFlags.add("REBOILER_TEMPERATURE_ABOVE_" + amineType + "_LIMIT");
    }
    if (overheadCO2Purity < 0.60) {
      degradationWarningFlags.add("HIGH_WATER_LOSS_LOW_OVERHEAD_CO2_PURITY");
    }
    if (leanLoading > 0.12) {
      degradationWarningFlags.add("LEAN_LOADING_ABOVE_TYPICAL_SPECIFICATION");
    }
    if (!electrolyteFlashSucceeded) {
      degradationWarningFlags.add("ELECTROLYTE_FLASH_NOT_AVAILABLE_CORRELATION_MODE");
    }
  }

  private double maxRecommendedReboilerTemperature() {
    return "MEA".equalsIgnoreCase(amineType) ? 122.0 : 130.0;
  }

  private double estimateLoading(SystemInterface system, double fallback) {
    double acidGas = getComponentMoles(system, "CO2") + getComponentMoles(system, "HCO3-") + getComponentMoles(system, "CO3--");
    double amine = getComponentMoles(system, amineType) + getComponentMoles(system, amineType + "+");
    return amine > 1.0e-12 ? acidGas / amine : fallback;
  }

  private double getComponentMoles(SystemInterface system, String componentName) {
    if (!system.hasComponent(componentName)) {
      return 0.0;
    }
    return system.getComponent(componentName).getNumberOfmoles();
  }

  private void removeComponentMoles(SystemInterface system, String componentName, double moles) {
    if (moles > 0.0 && system.hasComponent(componentName)) {
      system.addComponent(componentName, -Math.min(moles, getComponentMoles(system, componentName) * 0.999999));
    }
  }

  private static double limit(double value, double low, double high) {
    return Math.max(low, Math.min(high, value));
  }

  private static double toKw(double value, String unit) {
    if ("MW".equals(unit)) {
      return value * 1000.0;
    }
    return value;
  }

  private static double fromKw(double value, String unit) {
    if ("MW".equals(unit)) {
      return value / 1000.0;
    }
    return value;
  }

  private static double toBar(double value, String unit) {
    if ("Pa".equals(unit)) {
      return value / 1.0e5;
    }
    return value;
  }
}
