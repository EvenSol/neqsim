/*
 * Conductivity.java
 *
 * Created on 1. november 2000, 19:00
 */

package neqsim.physicalProperties.physicalPropertyMethods.liquidPhysicalProperties.diffusivity;
/**
 *
 * @author  Even Solbraa
 * @version
 */
public class CO2water extends Diffusivity{

    private static final long serialVersionUID = 1000;
    
    /** Creates new Conductivity */
    
    public CO2water() {
    }
    
    public CO2water(neqsim.physicalProperties.physicalPropertySystem.PhysicalPropertiesInterface liquidPhase) {
        super(liquidPhase);
    }
    
    // aqueous correlation
    public double calcBinaryDiffusionCoefficient(int i, int j, int method){
        binaryDiffusionCoeffisients[i][j] = 0.03389*Math.exp(-2213.7/liquidPhase.getPhase().getTemperature())*1e-4; // Tammi (1994) - Pcheco
        return binaryDiffusionCoeffisients[i][j];
    }
}
