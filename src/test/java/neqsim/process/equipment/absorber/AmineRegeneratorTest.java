package neqsim.process.equipment.absorber;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import neqsim.process.equipment.stream.Stream;
import neqsim.process.equipment.stream.StreamInterface;
import neqsim.thermo.system.SystemSrkEos;

/** Tests for amine solvent regeneration unit. */
public class AmineRegeneratorTest {
  @Test
  public void testMeaFlueGasRegeneration() {
    StreamInterface richMea = createAmineStream("MEA", 100.0, 40.0, 9.0, 1.0, 0.45);
    AmineRegenerator regenerator = new AmineRegenerator("MEA regenerator", richMea);
    regenerator.setAmineType("MEA");
    regenerator.setLeanLoading(0.08);
    regenerator.setStripperPressure(1.7, "bara");
    regenerator.setNumberOfStages(16);

    regenerator.run();

    assertNotNull(regenerator.getLeanAmineBottoms());
    assertNotNull(regenerator.getOverheadVapor());
    assertTrue(regenerator.getLeanLoading() < regenerator.getRichLoading());
    assertTrue(regenerator.getReboilerDuty("kW") > 0.0);
    assertTrue(regenerator.getCondenserDuty("kW") >= 0.0);
    assertTrue(regenerator.getOverheadCO2Purity() > 0.0);
    assertTrue(regenerator.getWaterMakeupRate("kg/hr") > 0.0);
  }

  @Test
  public void testMdeaNaturalGasSweeteningRegenerationAndAbsorberIntegration() {
    StreamInterface richMdea = createAmineStream("MDEA", 80.0, 55.0, 8.0, 1.0, 0.35);
    AmineRegenerator regenerator = new AmineRegenerator("MDEA regenerator", richMdea);
    regenerator.setAmineType("MDEA");
    regenerator.setCO2Recovery(0.85);
    regenerator.setPackingHeight(10.0);

    regenerator.run();

    RateBasedAbsorber absorber = new RateBasedAbsorber("integrated absorber");
    absorber.addSolventInStream(richMdea);
    regenerator.integrateWith(absorber);

    assertNotNull(absorber.getRegenerator());
    assertTrue(regenerator.getCO2Recovery() > 0.8);
    assertTrue(regenerator.getLeanLoading() < 0.10);
    assertTrue(regenerator.getSteamRate("kg/hr") > 0.0);
  }

  private StreamInterface createAmineStream(String amine, double amineMoles, double waterMoles, double temperatureC,
      double pressureBar, double loading) {
    SystemSrkEos system = new SystemSrkEos(273.15 + temperatureC, pressureBar);
    system.addComponent(amine, amineMoles);
    system.addComponent("water", waterMoles);
    system.addComponent("CO2", amineMoles * loading);
    system.setMixingRule("classic");
    Stream stream = new Stream("rich " + amine, system);
    stream.setFlowRate(1000.0, "kg/hr");
    return stream;
  }
}
