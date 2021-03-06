package neqsim.standards.util.example;

import neqsim.standards.StandardInterface;
import neqsim.standards.gasQuality.Standard_ISO6578;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicOperations.ThermodynamicOperations;

/*
 * PhaseEnvelope.java
 *
 * Created on 27. september 2001, 10:21
 */

/**
 *
 * @author esol
 * @version
 */
public class Test_ISO6578 {

    private static final long serialVersionUID = 1000;

    /**
     * Creates new PhaseEnvelope
     */
    public Test_ISO6578() {
    }

    public static void main(String args[]) {
        SystemInterface testSystem = new SystemSrkEos(273.15 - 160.0, 1.0);

        ThermodynamicOperations testOps = new ThermodynamicOperations(testSystem);
        testSystem.addComponent("nitrogen", 0.006538);
        testSystem.addComponent("methane", 0.91863);
        testSystem.addComponent("ethane", 0.058382);
        testSystem.addComponent("propane", 0.011993);
        // testSystem.addComponent("i-butane", 0.00);
        testSystem.addComponent("n-butane", 0.003255);
        testSystem.addComponent("i-pentane", 0.000657);
        testSystem.addComponent("n-pentane", 0.000545);
        
        

        testSystem.createDatabase(true);
        //testSystem.setMultiphaseWaxCheck(true);
        testSystem.setMixingRule(2);

        testSystem.init(0);
        StandardInterface standard = new Standard_ISO6578();//testSystem);
        standard.setThermoSystem(testSystem);
        standard.calculate();
        testSystem.display();

        System.out.println("corrfactor " + ((Standard_ISO6578) standard).getCorrFactor1());
      //  ((Standard_ISO6578) standard).useISO6578VolumeCorrectionFacotrs(false);

        standard.calculate();
        System.out.println("corrfactor " + ((Standard_ISO6578) standard).getCorrFactor1());
        testSystem.display();
    }
}
