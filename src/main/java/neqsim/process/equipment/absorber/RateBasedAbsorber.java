package neqsim.process.equipment.absorber;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import neqsim.process.equipment.stream.StreamInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermo.util.amines.AmineHeatOfAbsorption;

/**
 * Rate-based packed column absorber using the two-film model with mass transfer coefficients.
 *
 * <p>
 * This class implements a rigorous rate-based absorber/stripper model using the two-film theory with Onda (1968) or
 * Billet-Schultes (1999) correlations for mass transfer coefficients. It supports chemical enhancement factors for
 * reactive absorption (e.g., amine-CO2, amine-H2S).
 * </p>
 *
 * <p>
 * The column is divided into stages (segments), and for each stage the model:
 * </p>
 * <ol>
 * <li>Calculates gas-side and liquid-side mass transfer coefficients (kG, kL)</li>
 * <li>Applies enhancement factors for chemical reactions</li>
 * <li>Computes interfacial mass transfer rates using the overall driving force</li>
 * <li>Performs material and energy balances to update compositions and temperatures</li>
 * </ol>
 *
 * <h2>Mass Transfer Correlations</h2>
 * <p>
 * Supported mass transfer models via {@link MassTransferModel}:
 * </p>
 * <ul>
 * <li><b>ONDA_1968</b> — Onda, Takeuchi &amp; Okumoto, J. Chem. Eng. Japan, 1(1), 56-62 (1968)</li>
 * <li><b>BILLET_SCHULTES_1999</b> — Billet &amp; Schultes, Trans IChemE, 77(A), 498-504 (1999)</li>
 * </ul>
 *
 * <h2>Enhancement Factor Models</h2>
 * <ul>
 * <li><b>NONE</b> — Physical absorption only (E = 1.0)</li>
 * <li><b>HATTA_PSEUDO_FIRST_ORDER</b> — E = Ha/tanh(Ha) for pseudo-first-order reactions</li>
 * <li><b>VAN_KREVELEN_HOFTIJZER</b> — Generalized enhancement for fast reactions with finite reactant</li>
 * </ul>
 *
 * <h2>Usage Example</h2>
 *
 * <pre>
 * RateBasedAbsorber absorber = new RateBasedAbsorber("CO2 Absorber");
 * absorber.addGasInStream(sourGas);
 * absorber.addSolventInStream(leanAmine);
 * absorber.setNumberOfStages(20);
 * absorber.setPackedHeight(15.0);
 * absorber.setColumnDiameter(2.5);
 * absorber.setMassTransferModel(RateBasedAbsorber.MassTransferModel.ONDA_1968);
 * absorber.setEnhancementModel(RateBasedAbsorber.EnhancementModel.HATTA_PSEUDO_FIRST_ORDER);
 * absorber.setReactionRateConstant(5000.0); // 1/s for CO2-MDEA
 * absorber.run();
 * </pre>
 *
 * @author NeqSim
 * @version 1.0
 * @see SimpleAbsorber
 * @see SimpleTEGAbsorber
 */
public class RateBasedAbsorber extends SimpleAbsorber {
  /** Serialization version UID. */
  private static final long serialVersionUID = 1000;
  /** Logger object for class. */
  private static final Logger logger = LogManager.getLogger(RateBasedAbsorber.class);

  /**
   * Mass transfer correlation model.
   */
  public enum MassTransferModel {
    /**
     * Onda, Takeuchi &amp; Okumoto (1968). Standard correlation for random and structured packing. Applicable for
     * Reynolds numbers 0.04-500 (liquid), Schmidt numbers 0.1-10000.
     */
    ONDA_1968,

    /**
     * Billet &amp; Schultes (1999). More accurate for structured packing. Uses packing-specific constants CL and CV
     * from vendor data.
     */
    BILLET_SCHULTES_1999
  }

  /**
   * Enhancement factor model for chemical absorption.
   */
  public enum EnhancementModel {
    /** No enhancement — physical absorption only (E = 1.0). */
    NONE,

    /**
     * Pseudo-first-order approximation. E = Ha / tanh(Ha) where Ha = sqrt(k2 * D_A * C_B) / kL. Valid when Hatta number
     * is large and reactant is in excess.
     */
    HATTA_PSEUDO_FIRST_ORDER,

    /**
     * Van Krevelen-Hoftijzer model for fast reactions with finite reactant concentration. Accounts for depletion of
     * reactant B at the interface.
     */
    VAN_KREVELEN_HOFTIJZER
  }

  // ======================== Column geometry ========================
  /** Column internal diameter [m]. */
  private double columnDiameter = 1.0;

  /** Total packed height [m]. */
  private double packedHeight = 10.0;

  /** Packing specific surface area [m2/m3]. */
  private double packingSpecificArea = 250.0;

  /** Packing void fraction [-]. */
  private double packingVoidFraction = 0.95;

  /** Packing nominal size [m]. */
  private double packingNominalSize = 0.05;

  /** Critical surface tension of packing material [N/m]. */
  private double packingCriticalSurfaceTension = 0.075;

  // ======================== Model configuration ========================
  /** Mass transfer correlation. */
  private MassTransferModel massTransferModel = MassTransferModel.ONDA_1968;

  /** Enhancement factor model. */
  private EnhancementModel enhancementModel = EnhancementModel.NONE;

  /** Second-order reaction rate constant [m3/(mol·s)] for enhancement factor. */
  private double reactionRateConstant = 0.0;

  /** Stoichiometric ratio (moles of B consumed per mole of A absorbed). */
  private double stoichiometricRatio = 1.0;

  /** Billet-Schultes packing constant CL [-]. */
  private double billetCL = 1.0;

  /** Billet-Schultes packing constant CV [-]. */
  private double billetCV = 0.4;

  /** Explicit list of acid-gas components allowed to transfer. */
  private Set<String> transferableComponents = new LinkedHashSet<String>(Arrays.asList("CO2", "H2S"));

  /** Number of counter-current Picard iterations. */
  private int solverIterations = 6;

  // ======================== Operating parameters ========================
  /** Column operating pressure [Pa]. */
  private double operatingPressure = 101325.0;

  // ======================== Streams ========================
  /** Gas inlet stream (bottom of column). */
  private StreamInterface gasInStream;

  /** Solvent inlet stream (top of column). */
  private StreamInterface solventInStream;

  /** Gas outlet stream (top of column). */
  private StreamInterface gasOutStream;

  /** Solvent outlet stream (bottom of column). */
  private StreamInterface solventOutStream;

  // ======================== Results ========================
  /** Stage-by-stage results. */
  private List<StageResult> stageResults = new ArrayList<StageResult>();

  /** Overall gas-phase mass transfer coefficient times interfacial area [1/s]. */
  private double overallKGa = 0.0;

  /** Overall liquid-phase mass transfer coefficient times interfacial area [1/s]. */
  private double overallKLa = 0.0;

  /** Number of overall gas-phase transfer units [-]. */
  private double numberOfTransferUnits = 0.0;

  /** Height of an overall gas-phase transfer unit [m]. */
  private double heightOfTransferUnit = 0.0;

  /** Wetted specific area [m2/m3]. */
  private double wettedArea = 0.0;

  /** Column cross-sectional area [m2]. */
  private double columnArea = 0.0;

  /** Largest calculated temperature rise relative to lean solvent feed [K]. */
  private double temperatureBulge = 0.0;

  /**
   * Constructor for RateBasedAbsorber.
   *
   * @param name name of the absorber
   */
  public RateBasedAbsorber(String name) {
    super(name);
    setNumberOfStages(10);
  }

  /**
   * Set the gas inlet stream (enters at column bottom).
   *
   * @param stream gas inlet stream
   */
  public void addGasInStream(StreamInterface stream) {
    this.gasInStream = stream;
  }

  /**
   * Set the solvent inlet stream (enters at column top).
   *
   * @param stream solvent inlet stream
   */
  public void addSolventInStream(StreamInterface stream) {
    this.solventInStream = stream;
  }

  /**
   * Get the gas outlet stream (exits at column top).
   *
   * @return gas outlet stream
   */
  @Override
  public StreamInterface getGasOutStream() {
    if (gasOutStream == null && gasInStream != null) {
      gasOutStream = gasInStream.clone();
    }
    return gasOutStream;
  }

  /**
   * Get the solvent outlet stream (exits at column bottom).
   *
   * @return solvent outlet stream
   */
  public StreamInterface getSolventOutStream() {
    if (solventOutStream == null && solventInStream != null) {
      solventOutStream = solventInStream.clone();
    }
    return solventOutStream;
  }

  /**
   * Set the column internal diameter.
   *
   * @param diameter column diameter in meters
   */
  public void setColumnDiameter(double diameter) {
    this.columnDiameter = diameter;
  }

  /**
   * Get the column internal diameter.
   *
   * @return column diameter in meters
   */
  public double getColumnDiameter() {
    return columnDiameter;
  }

  /**
   * Set the total packed height.
   *
   * @param height packed height in meters
   */
  public void setPackedHeight(double height) {
    this.packedHeight = height;
  }

  /**
   * Get the total packed height.
   *
   * @return packed height in meters
   */
  public double getPackedHeight() {
    return packedHeight;
  }

  /**
   * Set the packing specific surface area.
   *
   * @param area specific surface area in m2/m3
   */
  public void setPackingSpecificArea(double area) {
    this.packingSpecificArea = area;
  }

  /**
   * Set the packing void fraction.
   *
   * @param voidFrac void fraction (0 to 1)
   */
  public void setPackingVoidFraction(double voidFrac) {
    this.packingVoidFraction = voidFrac;
  }

  /**
   * Set the packing nominal size.
   *
   * @param size nominal size in meters
   */
  public void setPackingNominalSize(double size) {
    this.packingNominalSize = size;
  }

  /**
   * Set the critical surface tension of packing material.
   *
   * @param sigma critical surface tension in N/m
   */
  public void setPackingCriticalSurfaceTension(double sigma) {
    this.packingCriticalSurfaceTension = sigma;
  }

  /**
   * Set the mass transfer correlation model.
   *
   * @param model mass transfer model to use
   */
  public void setMassTransferModel(MassTransferModel model) {
    this.massTransferModel = model;
  }

  /**
   * Get the mass transfer correlation model.
   *
   * @return current mass transfer model
   */
  public MassTransferModel getMassTransferModel() {
    return massTransferModel;
  }

  /**
   * Set the enhancement factor model.
   *
   * @param model enhancement model to use
   */
  public void setEnhancementModel(EnhancementModel model) {
    this.enhancementModel = model;
  }

  /**
   * Set the second-order reaction rate constant for enhancement calculations.
   *
   * @param k2 reaction rate constant in m3/(mol*s)
   */
  public void setReactionRateConstant(double k2) {
    this.reactionRateConstant = k2;
  }

  /**
   * Set the stoichiometric ratio for the absorption reaction.
   *
   * @param ratio moles of solvent reactant per mole of absorbed component
   */
  public void setStoichiometricRatio(double ratio) {
    this.stoichiometricRatio = ratio;
  }

  /**
   * Set Billet-Schultes packing constants.
   *
   * @param cl liquid-side packing constant CL
   * @param cv vapor-side packing constant CV
   */
  public void setBilletSchultesConstants(double cl, double cv) {
    this.billetCL = cl;
    this.billetCV = cv;
  }

  /**
   * Set the only components allowed to transfer between phases.
   *
   * @param componentNames component names, defaults should typically include CO2 and H2S
   */
  public void setTransferableComponents(String... componentNames) {
    transferableComponents.clear();
    if (componentNames != null) {
      for (String name : componentNames) {
        if (name != null && !name.trim().isEmpty()) {
          transferableComponents.add(name);
        }
      }
    }
  }

  /**
   * Get configured transferable components.
   *
   * @return copy of transferable component names
   */
  public Set<String> getTransferableComponents() {
    return new LinkedHashSet<String>(transferableComponents);
  }

  /**
   * Set number of counter-current solver iterations.
   *
   * @param iterations iterations, at least one
   */
  public void setSolverIterations(int iterations) {
    solverIterations = Math.max(1, iterations);
  }

  /**
   * Get the number of overall gas-phase transfer units.
   *
   * @return NTU value
   */
  public double getNumberOfTransferUnits() {
    return numberOfTransferUnits;
  }

  /**
   * Get the height of an overall gas-phase transfer unit.
   *
   * @return HTU in meters
   */
  public double getHeightOfTransferUnit() {
    return heightOfTransferUnit;
  }

  /**
   * Get the overall gas-side mass transfer coefficient times interfacial area.
   *
   * @return KGa in 1/s
   */
  public double getOverallKGa() {
    return overallKGa;
  }

  /**
   * Get the overall liquid-side mass transfer coefficient times interfacial area.
   *
   * @return KLa in 1/s
   */
  public double getOverallKLa() {
    return overallKLa;
  }

  /**
   * Get the wetted specific area of the packing.
   *
   * @return wetted area in m2/m3
   */
  @Override
  public double getWettedArea() {
    return wettedArea;
  }

  /**
   * Get the column temperature bulge relative to lean solvent inlet.
   *
   * @return temperature bulge in kelvin
   */
  public double getTemperatureBulge() {
    return temperatureBulge;
  }

  /**
   * Get the stage-by-stage calculation results.
   *
   * @return list of stage results
   */
  public List<StageResult> getStageResults() {
    return stageResults;
  }

  /** {@inheritDoc} */
  @Override
  public void run(UUID id) {
    if (gasInStream == null || solventInStream == null) {
      logger.error("Both gas and solvent inlet streams must be set");
      return;
    }

    columnArea = Math.PI / 4.0 * columnDiameter * columnDiameter;
    double stageHeight = packedHeight / getNumberOfStages();
    stageResults.clear();
    temperatureBulge = 0.0;

    SystemInterface[] gasProfile = new SystemInterface[getNumberOfStages() + 1];
    SystemInterface[] liquidProfile = new SystemInterface[getNumberOfStages() + 1];
    gasProfile[0] = initializedClone(gasInStream.getThermoSystem());
    liquidProfile[getNumberOfStages()] = initializedClone(solventInStream.getThermoSystem());
    operatingPressure = gasProfile[0].getPressure("Pa");

    for (int i = 1; i <= getNumberOfStages(); i++) {
      gasProfile[i] = initializedClone(gasProfile[0]);
    }
    for (int i = 0; i < getNumberOfStages(); i++) {
      liquidProfile[i] = initializedClone(liquidProfile[getNumberOfStages()]);
    }

    for (int iteration = 0; iteration < solverIterations; iteration++) {
      SystemInterface[] nextGas = new SystemInterface[getNumberOfStages() + 1];
      SystemInterface[] nextLiquid = new SystemInterface[getNumberOfStages() + 1];
      nextGas[0] = initializedClone(gasProfile[0]);
      nextLiquid[getNumberOfStages()] = initializedClone(liquidProfile[getNumberOfStages()]);
      List<StageResult> iterationResults = new ArrayList<StageResult>();

      for (int stage = 0; stage < getNumberOfStages(); stage++) {
        SystemInterface gasIn = initializedClone(nextGas[stage]);
        SystemInterface liquidIn = initializedClone(liquidProfile[stage + 1]);
        StageResult result = calculateStage(gasIn, liquidIn, stageHeight, stage);
        nextGas[stage + 1] = gasIn;
        nextLiquid[stage] = liquidIn;
        iterationResults.add(result);
      }

      gasProfile = nextGas;
      liquidProfile = nextLiquid;
      stageResults = iterationResults;
    }

    double avgKGa = 0.0;
    double avgKLa = 0.0;
    double avgWettedArea = 0.0;
    double leanTemperature = solventInStream.getThermoSystem().getTemperature();
    for (StageResult result : stageResults) {
      avgKGa += result.kGa;
      avgKLa += result.kLa;
      avgWettedArea += result.wettedArea;
      temperatureBulge = Math.max(temperatureBulge, result.temperature - leanTemperature);
    }

    overallKGa = avgKGa / getNumberOfStages();
    overallKLa = avgKLa / getNumberOfStages();
    wettedArea = avgWettedArea / getNumberOfStages();

    if (overallKGa > 0.0) {
      double gasVelocity = gasInStream.getFlowRate("kg/hr") / 3600.0 / gasProfile[0].getDensity("kg/m3") / columnArea;
      heightOfTransferUnit = gasVelocity / overallKGa;
      if (heightOfTransferUnit > 0.0) {
        numberOfTransferUnits = packedHeight / heightOfTransferUnit;
      }
    }

    gasOutStream = gasInStream.clone();
    gasOutStream.setThermoSystem(gasProfile[getNumberOfStages()]);
    gasOutStream.setCalculationIdentifier(id);

    solventOutStream = solventInStream.clone();
    solventOutStream.setThermoSystem(liquidProfile[0]);
    solventOutStream.setCalculationIdentifier(id);

    setCalculationIdentifier(id);
  }

  private SystemInterface initializedClone(SystemInterface source) {
    SystemInterface clone = source.clone();
    ThermodynamicOperations ops = new ThermodynamicOperations(clone);
    try {
      ops.TPflash();
    } catch (Exception ex) {
      logger.warn("Rate-based absorber flash warning: {}", ex.getMessage());
    }
    clone.initProperties();
    return clone;
  }

  /**
   * Calculate mass transfer for a single stage.
   *
   * @param gasPhase gas-phase system (modified in place)
   * @param liquidPhase liquid-phase system (modified in place)
   * @param stageHeight height of this stage segment [m]
   * @param stageIndex stage number (0-based)
   * @return stage calculation result
   */
  private StageResult calculateStage(SystemInterface gasPhase, SystemInterface liquidPhase, double stageHeight,
      int stageIndex) {

    StageResult result = new StageResult();
    result.stageNumber = stageIndex + 1;
    result.temperature = gasPhase.getTemperature();
    result.pressure = gasPhase.getPressure("Pa");

    // Get phase properties
    double rhoG = gasPhase.getDensity("kg/m3");
    double rhoL = liquidPhase.getDensity("kg/m3");
    double muG = gasPhase.getPhase(0).getViscosity("kg/msec");
    double muL = liquidPhase.getPhase(0).getViscosity("kg/msec");
    double sigmaL = gasPhase.getInterphaseProperties().getSurfaceTension(0, 1);
    if (sigmaL <= 0.0) {
      sigmaL = 0.03; // Default if not available
    }

    // Superficial velocities
    double gasFlow = gasInStream.getFlowRate("kg/hr") / 3600.0;
    double liquidFlow = solventInStream.getFlowRate("kg/hr") / 3600.0;
    double uG = gasFlow / (rhoG * columnArea);
    double uL = liquidFlow / (rhoL * columnArea);

    double dG = componentDiffusivity(gasPhase, "CO2", 1.5e-5, true);
    double dL = componentDiffusivity(liquidPhase, "CO2", 1.5e-9, false);

    // Calculate mass transfer coefficients based on selected model
    double kG = 0.0;
    double kL = 0.0;
    double aw = 0.0;

    if (massTransferModel == MassTransferModel.ONDA_1968) {
      double[] ondaResult = calculateOndaMassTransfer(uG, uL, rhoG, rhoL, muG, muL, sigmaL, dG, dL);
      kG = ondaResult[0];
      kL = ondaResult[1];
      aw = ondaResult[2];
    } else {
      double[] billetResult = calculateBilletSchultesMassTransfer(uG, uL, rhoG, rhoL, muG, muL, sigmaL, dG, dL);
      kG = billetResult[0];
      kL = billetResult[1];
      aw = billetResult[2];
    }

    // Apply enhancement factor for chemical absorption
    double enhancementFactor = calculateEnhancementFactor(kL, liquidPhase, rhoL, muL, dL);

    double kLEnhanced = kL * enhancementFactor;
    result.enhancementFactor = enhancementFactor;
    result.kGa = kG * aw;
    result.kLa = kLEnhanced * aw;
    result.wettedArea = aw;

    // Calculate mass transfer for each transferring component
    double totalMolesTransferred = 0.0;
    int nComp = gasPhase.getNumberOfComponents();
    for (int i = 0; i < nComp; i++) {
      String compName = gasPhase.getComponent(i).getName();
      // Only explicitly configured acid gases are transferable.
      if (!transferableComponents.contains(compName) || !liquidPhase.hasComponent(compName)) {
        continue;
      }

      // Mole fractions
      double yBulk = gasPhase.getPhase(0).getComponent(i).getx();
      double xBulk = 0.0;
      int liqIndex = liquidPhase.getPhase(0).getComponent(compName).getComponentNumber();
      if (liqIndex >= 0) {
        xBulk = liquidPhase.getPhase(0).getComponent(liqIndex).getx();
      }

      // Equilibrium: y* = K * x (VLE K-value)
      double Ki = gasPhase.getPhase(0).getComponent(i).getFugacityCoefficient()
          / liquidPhase.getPhase(0).getComponent(liqIndex).getFugacityCoefficient();
      if (Ki <= 0.0) {
        Ki = 1.0;
      }
      double yEquil = Ki * xBulk;

      // Overall driving force (gas-side)
      double drivingForce = yBulk - yEquil;
      result.approachToEquilibrium = drivingForce;
      if (Math.abs(drivingForce) < 1e-15) {
        continue;
      }

      // Overall gas-phase mass transfer coefficient: 1/KOG = 1/kG + m/kL_enhanced
      // where m = slope of equilibrium line (Henry's law constant or K-value)
      double mSlope = Ki;
      double koG = 0.0;
      if (kG > 0 && kLEnhanced > 0) {
        koG = 1.0 / (1.0 / kG + mSlope / kLEnhanced);
      }

      // Molar flux [mol/(m2·s)]
      double totalMolarDensity = gasPhase.getPressure("Pa") / (8.314 * gasPhase.getTemperature());
      double flux = koG * totalMolarDensity * drivingForce;

      // Total moles transferred in this stage
      double interfacialArea = aw * columnArea * stageHeight;
      double molesTransferred = flux * interfacialArea;

      // Update compositions: remove from gas, add to liquid
      if (molesTransferred > 0) {
        // Positive = transfer from gas to liquid (absorption)
        double maxRemovable = gasPhase.getComponent(i).getNumberOfmoles() * 0.9;
        double actualTransfer = Math.min(molesTransferred, maxRemovable);
        if (actualTransfer > 1e-20) {
          gasPhase.addComponent(i, -actualTransfer);
          liquidPhase.addComponent(liqIndex, actualTransfer);
          totalMolesTransferred += actualTransfer;
          result.heatReleased += calculateHeatReleased(compName, actualTransfer, liquidPhase);
        }
      }
    }

    result.molesTransferred = totalMolesTransferred;
    applyTemperatureRise(gasPhase, liquidPhase, result.heatReleased);

    // Re-flash both phases after composition update
    ThermodynamicOperations gasOps = new ThermodynamicOperations(gasPhase);
    try {
      gasOps.TPflash();
    } catch (Exception ex) {
      logger.warn("Stage {} gas flash warning: {}", stageIndex, ex.getMessage());
    }
    gasPhase.initProperties();

    ThermodynamicOperations liqOps = new ThermodynamicOperations(liquidPhase);
    try {
      liqOps.TPflash();
    } catch (Exception ex) {
      logger.warn("Stage {} liquid flash warning: {}", stageIndex, ex.getMessage());
    }
    liquidPhase.initProperties();

    result.temperature = 0.5 * (gasPhase.getTemperature() + liquidPhase.getTemperature());
    result.loading = calculateLoading(liquidPhase);
    if (gasPhase.hasComponent("CO2")) {
      result.co2MoleFraction = gasPhase.getPhase(0).getComponent("CO2").getx();
    }

    return result;
  }

  /**
   * Calculate mass transfer coefficients using the Onda (1968) correlation.
   *
   * <p>
   * Gas-side: kG = C1 * (a_p * D_G) * (Re_G)^0.7 * (Sc_G)^(1/3) * (a_p * d_p)^(-2.0)
   * </p>
   * <p>
   * Liquid-side: kL = 0.0051 * (Re_L)^(2/3) * (Sc_L)^(-0.5) * (a_p * d_p)^0.4 * (mu_L * g/rho_L)^(1/3)
   * </p>
   * <p>
   * Wetted area: a_w/a_p = 1 - exp(-1.45 * (sigma_c/sigma_L)^0.75 * Re_L^0.1 * Fr_L^(-0.05) * We_L^0.2)
   * </p>
   *
   * @param uG gas superficial velocity [m/s]
   * @param uL liquid superficial velocity [m/s]
   * @param rhoG gas density [kg/m3]
   * @param rhoL liquid density [kg/m3]
   * @param muG gas viscosity [Pa.s]
   * @param muL liquid viscosity [Pa.s]
   * @param sigmaL liquid surface tension [N/m]
   * @return array [kG, kL, aw] where kG in [mol/(m2.s.Pa)], kL in [m/s], aw in [m2/m3]
   */
  private double[] calculateOndaMassTransfer(double uG, double uL, double rhoG, double rhoL, double muG, double muL,
      double sigmaL, double dG, double dL) {
    double g = 9.81;
    double ap = packingSpecificArea;
    double dp = packingNominalSize;
    double sigmaC = packingCriticalSurfaceTension;

    // Reynolds numbers
    double reG = uG * rhoG / (muG * ap);
    double reL = uL * rhoL / (muL * ap);
    reL = Math.max(reL, 1e-6);

    // Schmidt numbers
    double scG = muG / (rhoG * dG);
    double scL = muL / (rhoL * dL);

    // Froude and Weber numbers for wetted area
    double frL = uL * uL * ap / g;
    double weL = uL * uL * rhoL / (sigmaL * ap);

    // Wetted area ratio (Onda 1968)
    double sigmaRatio = sigmaC / sigmaL;
    sigmaRatio = Math.min(sigmaRatio, 1.0); // Cap at 1.0 per Onda
    double awRatio = 1.0 - Math.exp(-1.45 * Math.pow(sigmaRatio, 0.75) * Math.pow(Math.max(reL, 0.01), 0.1)
        * Math.pow(Math.max(frL, 1e-6), -0.05) * Math.pow(Math.max(weL, 1e-6), 0.2));
    awRatio = Math.max(awRatio, 0.1);
    awRatio = Math.min(awRatio, 1.0);
    double aw = awRatio * ap;

    // Gas-phase mass transfer coefficient (Onda 1968)
    // kG * R * T / (a_p * D_G) = C1 * Re_G^0.7 * Sc_G^(1/3) * (a_p * d_p)^(-2.0)
    double c1 = 5.23; // Onda constant for kG
    double kGDimensionless = c1 * Math.pow(Math.max(reG, 0.01), 0.7) * Math.pow(scG, 1.0 / 3.0)
        * Math.pow(ap * dp, -2.0);
    double kG = kGDimensionless * ap * dG; // [m/s] gas-side MTC

    // Liquid-phase mass transfer coefficient (Onda 1968)
    // kL * (rho_L / (mu_L * g))^(1/3) = 0.0051 * (Re_L)^(2/3) * (Sc_L)^(-0.5) * (a_p*d_p)^0.4
    double kLDimensionless = 0.0051 * Math.pow(Math.max(reL, 0.01), 2.0 / 3.0) * Math.pow(scL, -0.5)
        * Math.pow(ap * dp, 0.4);
    double kL = kLDimensionless * Math.pow(muL * g / rhoL, 1.0 / 3.0); // [m/s]

    return new double[] { kG, kL, aw };
  }

  /**
   * Calculate mass transfer coefficients using the Billet-Schultes (1999) correlation.
   *
   * @param uG gas superficial velocity [m/s]
   * @param uL liquid superficial velocity [m/s]
   * @param rhoG gas density [kg/m3]
   * @param rhoL liquid density [kg/m3]
   * @param muG gas viscosity [Pa.s]
   * @param muL liquid viscosity [Pa.s]
   * @param sigmaL liquid surface tension [N/m]
   * @return array [kG, kL, aw] where kG in [m/s], kL in [m/s], aw in [m2/m3]
   */
  private double[] calculateBilletSchultesMassTransfer(double uG, double uL, double rhoG, double rhoL, double muG,
      double muL, double sigmaL, double dG, double dL) {
    double g = 9.81;
    double ap = packingSpecificArea;
    double eps = packingVoidFraction;

    // Effective velocity
    double uGe = uG / (eps * (1 - calculateLiquidHoldup(uL, rhoL, muL)));
    double uLe = uL / (eps * calculateLiquidHoldup(uL, rhoL, muL));

    // Hydraulic diameter
    double dh = 4.0 * eps / ap;

    // Reynolds numbers
    double reG = uGe * dh * rhoG / muG;
    double reL = uLe * dh * rhoL / muL;
    reL = Math.max(reL, 1e-6);

    // Schmidt numbers
    double scG = muG / (rhoG * dG);
    double scL = muL / (rhoL * dL);

    // Billet-Schultes gas-side MTC
    double kG = billetCV / Math.pow(dh, 0.5) * Math.sqrt(dG) * Math.pow(Math.max(reG, 0.01), 0.75)
        * Math.pow(scG, 1.0 / 3.0);

    // Billet-Schultes liquid-side MTC
    double kL = billetCL * Math.pow(12.0, 1.0 / 6.0) * Math.pow(dL / dh, 0.5) * Math.pow(uLe, 1.0 / 3.0)
        * Math.pow(g, 1.0 / 6.0);

    // Wetted area (simplified — same as Onda for random, or vendor data)
    double aw = 0.85 * ap; // Default fraction for structured packing

    return new double[] { kG, kL, aw };
  }

  /**
   * Calculate operating liquid holdup (simplified Stichlmair model).
   *
   * @param uL liquid superficial velocity [m/s]
   * @param rhoL liquid density [kg/m3]
   * @param muL liquid viscosity [Pa.s]
   * @return liquid holdup fraction [-]
   */
  private double calculateLiquidHoldup(double uL, double rhoL, double muL) {
    double g = 9.81;
    double ap = packingSpecificArea;
    double reL = uL * rhoL / (muL * ap);
    reL = Math.max(reL, 1e-6);
    double hL = 0.555 * Math.pow(reL, 1.0 / 3.0) * Math.pow(ap, 2.0 / 3.0) / Math.pow(rhoL * g / muL, 1.0 / 3.0);
    hL = Math.max(hL, 0.01);
    hL = Math.min(hL, 0.3);
    return hL;
  }

  /**
   * Calculate enhancement factor for chemical absorption.
   *
   * @param kL physical liquid-side mass transfer coefficient [m/s]
   * @param liquidPhase liquid phase system
   * @param rhoL liquid density [kg/m3]
   * @param muL liquid viscosity [Pa.s]
   * @return enhancement factor E [-]
   */
  private double calculateEnhancementFactor(double kL, SystemInterface liquidPhase, double rhoL, double muL, double dL) {
    if (enhancementModel == EnhancementModel.NONE || reactionRateConstant <= 0.0) {
      return 1.0;
    }

    double dB = componentDiffusivity(liquidPhase, dominantAmineName(liquidPhase), 1.0e-9, false);

    if (enhancementModel == EnhancementModel.HATTA_PSEUDO_FIRST_ORDER) {
      // Pseudo-first-order: Ha = sqrt(k2 * C_B * D_A) / kL
      // E = Ha / tanh(Ha)
      double cB = calculateReactantConcentration(liquidPhase);
      double ha = Math.sqrt(reactionRateConstant * amineKineticMultiplier(liquidPhase) * cB * dL) / Math.max(kL, 1e-10);
      if (ha < 0.3) {
        return 1.0; // Slow reaction regime
      }
      return ha / Math.tanh(ha);
    }

    if (enhancementModel == EnhancementModel.VAN_KREVELEN_HOFTIJZER) {
      // Van Krevelen-Hoftijzer for finite reactant
      double cB = calculateReactantConcentration(liquidPhase);
      double ha = Math.sqrt(reactionRateConstant * amineKineticMultiplier(liquidPhase) * cB * dL) / Math.max(kL, 1e-10);
      // E_inf = 1 + (D_B/D_A) * (C_B / (nu * C_Ai))
      double eInf = 1.0 + (dB / dL) * cB / (stoichiometricRatio * 1.0);
      eInf = Math.max(eInf, 1.0);

      // Approximate: E ≈ min(E from Ha, E_inf)
      double eHa = ha / Math.tanh(ha);
      return Math.min(eHa, eInf);
    }

    return 1.0;
  }


  private double componentDiffusivity(SystemInterface system, String component, double fallback, boolean gasPhase) {
    if (component == null || !system.hasComponent(component)) {
      return fallback;
    }
    try {
      int componentIndex = system.getPhase(0).getComponent(component).getComponentNumber();
      double value = system.getPhase(0).getPhysicalProperties().getEffectiveDiffusionCoefficient(componentIndex);
      if (value > 0.0 && Double.isFinite(value)) {
        return value;
      }
    } catch (RuntimeException ex) {
      logger.debug("Using fallback {} diffusivity for {}: {}", gasPhase ? "gas" : "liquid", component,
          ex.getMessage());
    }
    return fallback;
  }

  private double calculateReactantConcentration(SystemInterface liquidPhase) {
    double amineMoles = 0.0;
    String[] amines = new String[] { "MEA", "DEA", "MDEA", "PZ", "piperazine", "aMDEA" };
    for (String amine : amines) {
      if (liquidPhase.hasComponent(amine)) {
        amineMoles += liquidPhase.getComponent(amine).getNumberOfmoles();
      }
    }
    double volume = Math.max(liquidPhase.getVolume("m3"), 1.0e-12);
    return Math.max(amineMoles / volume, 1.0e-9);
  }

  private String dominantAmineName(SystemInterface liquidPhase) {
    String dominant = "water";
    double maxMoles = 0.0;
    String[] amines = new String[] { "MEA", "DEA", "MDEA", "PZ", "piperazine", "aMDEA" };
    for (String amine : amines) {
      if (liquidPhase.hasComponent(amine)) {
        double moles = liquidPhase.getComponent(amine).getNumberOfmoles();
        if (moles > maxMoles) {
          maxMoles = moles;
          dominant = amine;
        }
      }
    }
    return dominant;
  }

  private double amineKineticMultiplier(SystemInterface liquidPhase) {
    if (liquidPhase.hasComponent("MEA")) {
      return 8.0;
    }
    if (liquidPhase.hasComponent("DEA")) {
      return 4.0;
    }
    if (liquidPhase.hasComponent("PZ") || liquidPhase.hasComponent("piperazine") || liquidPhase.hasComponent("aMDEA")) {
      return 12.0;
    }
    if (liquidPhase.hasComponent("MDEA")) {
      return 1.0;
    }
    return 0.0;
  }

  private double calculateHeatReleased(String componentName, double molesTransferred, SystemInterface liquidPhase) {
    if (molesTransferred <= 0.0) {
      return 0.0;
    }
    if ("CO2".equals(componentName) || "H2S".equals(componentName)) {
      AmineHeatOfAbsorption heat = new AmineHeatOfAbsorption(resolveAmineType(liquidPhase), calculateAmineMassFraction(liquidPhase),
          calculateLoading(liquidPhase), liquidPhase.getTemperature());
      double heatKjPerMol = "H2S".equals(componentName) ? heat.calcHeatOfAbsorptionH2S() : heat.calcHeatOfAbsorptionCO2();
      return Math.abs(heatKjPerMol) * 1000.0 * molesTransferred;
    }
    return 0.0;
  }

  private AmineHeatOfAbsorption.AmineType resolveAmineType(SystemInterface liquidPhase) {
    if (liquidPhase.hasComponent("MEA")) {
      return AmineHeatOfAbsorption.AmineType.MEA;
    }
    if (liquidPhase.hasComponent("DEA")) {
      return AmineHeatOfAbsorption.AmineType.DEA;
    }
    if (liquidPhase.hasComponent("PZ") || liquidPhase.hasComponent("piperazine") || liquidPhase.hasComponent("aMDEA")) {
      return AmineHeatOfAbsorption.AmineType.AMDEA;
    }
    return AmineHeatOfAbsorption.AmineType.MDEA;
  }

  private double calculateAmineMassFraction(SystemInterface liquidPhase) {
    double amineMass = 0.0;
    String[] amines = new String[] { "MEA", "DEA", "MDEA", "PZ", "piperazine", "aMDEA" };
    for (String amine : amines) {
      if (liquidPhase.hasComponent(amine)) {
        amineMass += liquidPhase.getComponent(amine).getNumberOfmoles() * liquidPhase.getComponent(amine).getMolarMass();
      }
    }
    double totalMass = Math.max(liquidPhase.getMass("kg"), 1.0e-12);
    return Math.max(0.0, Math.min(0.9, amineMass / totalMass));
  }

  private double calculateLoading(SystemInterface liquidPhase) {
    double acidMoles = 0.0;
    if (liquidPhase.hasComponent("CO2")) {
      acidMoles += liquidPhase.getComponent("CO2").getNumberOfmoles();
    }
    double amineMoles = Math.max(calculateReactantConcentration(liquidPhase) * liquidPhase.getVolume("m3"), 1.0e-12);
    return acidMoles / amineMoles;
  }

  private void applyTemperatureRise(SystemInterface gasPhase, SystemInterface liquidPhase, double heatReleased) {
    if (heatReleased <= 0.0) {
      return;
    }
    double cp = Math.max(Math.abs(gasPhase.getCp()) + Math.abs(liquidPhase.getCp()), 1.0e3);
    double deltaT = Math.min(25.0, heatReleased / cp);
    gasPhase.setTemperature(gasPhase.getTemperature() + deltaT);
    liquidPhase.setTemperature(liquidPhase.getTemperature() + deltaT);
  }

  /**
   * Stage-by-stage calculation result container.
   */
  public static class StageResult implements Serializable {
    /** Serialization version UID. */
    private static final long serialVersionUID = 1000;

    /** Stage number (1-based). */
    public int stageNumber;

    /** Stage temperature [K]. */
    public double temperature;

    /** Stage pressure [Pa]. */
    public double pressure;

    /** Gas-side KGa [1/s]. */
    public double kGa;

    /** Liquid-side KLa [1/s]. */
    public double kLa;

    /** Wetted area [m2/m3]. */
    public double wettedArea;

    /** Enhancement factor [-]. */
    public double enhancementFactor = 1.0;

    /** Total moles transferred in this stage [mol/s]. */
    public double molesTransferred;

    /** Heat released in this stage [W]. */
    public double heatReleased;

    /** CO2 mole fraction in exiting gas [-]. */
    public double co2MoleFraction;

    /** CO2 loading in exiting liquid [mol/mol amine]. */
    public double loading;

    /** Approach to equilibrium y-y* [-]. */
    public double approachToEquilibrium;

    /**
     * Get the stage number.
     *
     * @return stage number, one based
     */
    public int getStageNumber() {
      return stageNumber;
    }

    /**
     * Get the stage temperature.
     *
     * @return stage temperature in kelvin
     */
    public double getTemperature() {
      return temperature;
    }

    /**
     * Get the stage pressure.
     *
     * @return stage pressure in pascal
     */
    public double getPressure() {
      return pressure;
    }

    /**
     * Get the gas-side volumetric mass-transfer coefficient.
     *
     * @return gas-side KGa in 1/s
     */
    public double getKGa() {
      return kGa;
    }

    /**
     * Get the liquid-side volumetric mass-transfer coefficient.
     *
     * @return liquid-side KLa in 1/s
     */
    public double getKLa() {
      return kLa;
    }

    /**
     * Get the wetted packing area.
     *
     * @return wetted area in m2/m3
     */
    public double getWettedArea() {
      return wettedArea;
    }

    /**
     * Get the enhancement factor.
     *
     * @return reaction enhancement factor
     */
    public double getEnhancementFactor() {
      return enhancementFactor;
    }

    /**
     * Get the total moles transferred in this stage.
     *
     * @return moles transferred in mol/s
     */
    public double getMolesTransferred() {
      return molesTransferred;
    }

    public double getHeatReleased() {
      return heatReleased;
    }

    public double getCo2MoleFraction() {
      return co2MoleFraction;
    }

    public double getLoading() {
      return loading;
    }

    public double getApproachToEquilibrium() {
      return approachToEquilibrium;
    }
  }
}
