package neqsim.thermo.util.referenceequations;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/** Test class for {@link IAPWSIF97}. */
public class IAPWSIF97Test {

  @Test
  public void testSampleValues() {
    double h = IAPWSIF97.h_pT(3.0, 573.15);
    double v = IAPWSIF97.v_pT(3.0, 573.15);
    assertEquals(2994.3, h, 1.0);
    assertEquals(0.081, v, 0.05);

    double psat = IAPWSIF97.psat_T(373.15);
    assertEquals(0.1014, psat, 1e-3);

    double Tsat = IAPWSIF97.Tsat_p(1.0);
    assertEquals(453.0, Tsat, 1.0);

    assertTrue(!Double.isNaN(h));
  }
}
