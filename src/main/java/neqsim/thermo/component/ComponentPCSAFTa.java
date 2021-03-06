  /*
 * System_SRK_EOS.java
 *
 * Created on 8. april 2000, 23:14
 */
package neqsim.thermo.component;

import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseInterface;

/**
 *
 * @author  Even Solbraa
 * @version
 */
public class ComponentPCSAFTa extends ComponentPCSAFT implements ComponentCPAInterface {

    private static final long serialVersionUID = 1000;

    /** Creates new System_SRK_EOS
     * Ev liten fil ja.
     */
    int cpaon = 1;
    private double[][] xsitedni = new double[0][0];
    double[] xsite = new double[0];
    double[] xsiteOld = new double[0];
    double[] xsitedV = new double[0];
    double[] xsitedT = new double[0];

    public ComponentPCSAFTa() {
    }

    public ComponentPCSAFTa(double moles) {
        super(moles);
    }

    public ComponentPCSAFTa(String component_name, double moles, double molesInPhase, int compnumber) {
        super(component_name, moles, molesInPhase, compnumber);
        xsite = new double[numberOfAssociationSites];
        xsitedV = new double[numberOfAssociationSites];
        xsiteOld = new double[numberOfAssociationSites];
        xsitedT = new double[numberOfAssociationSites];

        if (numberOfAssociationSites != 0 && cpaon == 1) {
            for (int j = 0; j < getNumberOfAssociationSites(); j++) {
                setXsite(j, 1.0);
                setXsiteOld(j, 1.0);
                setXsitedV(j, 0.0);
                setXsitedT(j, 0.0);
            }
        }
    }

    public ComponentPCSAFTa(int number, double TC, double PC, double M, double a, double moles) {
        super(number, TC, PC, M, a, moles);
        xsite = new double[numberOfAssociationSites];
        xsitedV = new double[numberOfAssociationSites];
        xsiteOld = new double[numberOfAssociationSites];
        xsitedT = new double[numberOfAssociationSites];

        if (numberOfAssociationSites != 0 && cpaon == 1) {
            for (int j = 0; j < getNumberOfAssociationSites(); j++) {
                setXsite(j, 1.0);
                setXsiteOld(j, 1.0);
                setXsitedV(j, 0.0);
                setXsitedT(j, 0.0);
            }
        }
    }

    public Object clone() {

        ComponentPCSAFTa clonedComponent = null;
        try {
            clonedComponent = (ComponentPCSAFTa) super.clone();
        } catch (Exception e) {
            logger.error("Cloning failed.", e);
        }

        clonedComponent.xsite = xsite.clone();
        System.arraycopy(this.xsite, 0, clonedComponent.xsite, 0, xsite.length);
        clonedComponent.xsiteOld = xsiteOld.clone();
        System.arraycopy(this.xsiteOld, 0, clonedComponent.xsiteOld, 0, xsiteOld.length);
        clonedComponent.xsitedV = xsitedV.clone();
        System.arraycopy(this.xsitedV, 0, clonedComponent.xsitedV, 0, xsitedV.length);
        clonedComponent.xsitedT = xsitedT.clone();
        System.arraycopy(this.xsitedT, 0, clonedComponent.xsitedT, 0, xsitedT.length);
        return clonedComponent;
    }

    public void init(double temperature, double pressure, double totalNumberOfMoles, double beta, int type) {
        super.init(temperature, pressure, totalNumberOfMoles, beta, type);
    }

    public double dFdN(PhaseInterface phase, int numberOfComponents, double temperature, double pressure) {
        double Fsup = super.dFdN(phase, numberOfComponents, temperature, pressure);
        double Fcpa = 0.0;
        //if(phase.getPhaseType()==1) cpaon=0;
        Fcpa = dFCPAdN(phase, numberOfComponents, temperature, pressure);
        //System.out.println("Fsup " + Fsup + "  fcpa " + Fcpa);
        return Fsup + cpaon * Fcpa;
    }

    public double dFdNdT(PhaseInterface phase, int numberOfComponents, double temperature, double pressure) {
        return super.dFdNdT(phase, numberOfComponents, temperature, pressure);
    }

    public double dFdNdV(PhaseInterface phase, int numberOfComponents, double temperature, double pressure) {
        return super.dFdNdV(phase, numberOfComponents, temperature, pressure);
    }

    public double dFdNdN(int j, PhaseInterface phase, int numberOfComponents, double temperature, double pressure) {
        return super.dFdNdN(j, phase, numberOfComponents, temperature, pressure);
    }

    public double dFCPAdN(PhaseInterface phase, int numberOfComponents, double temperature, double pressure) {
        double xi = 0.0;
        for (int i = 0; i < numberOfAssociationSites; i++) {
            xi += Math.log(xsite[i]);
        }
        return (xi - ((PhaseCPAInterface) phase).getHcpatot() / 2.0 * dlogghsSAFTdi);//calc_lngi(phase));
    }

    public double dFCPAdNdV(PhaseInterface phase, int numberOfComponents, double temperature, double pressure) {
        double xi = dFCPAdNdXidXdV(phase);
        double xi2 = -((PhaseCPAInterface) phase).getHcpatot() / 2.0 * calc_lngidV(phase);
        return xi + xi2;
    }

    public double calc_lngidV(PhaseInterface phase) {
        return 0;
    }

    public double dFCPAdVdXi(int site, PhaseInterface phase) {
        return 1.0 / (2.0 * phase.getTotalVolume()) * (1.0 - phase.getTotalVolume() * ((PhaseCPAInterface) phase).getGcpav()) * getNumberOfMolesInPhase();
    }

    public double dFCPAdNdXi(int site, PhaseInterface phase) {
        double xi = 1.0 / xsite[site];

        //  return xi - tempp;
        return xi + getNumberOfMolesInPhase() / 2.0 * calc_lngi(phase);
    }

    public double dFCPAdXidXj(int sitei, int sitej, int compj, PhaseInterface phase) {
        return 0.0;
    }

    public double dFCPAdXi(int site, PhaseInterface phase) {
        return 0.0;
    }

    public double calc_lngi(PhaseInterface phase) {
        return 0;
    }

    public double dFCPAdNdXidXdV(PhaseInterface phase) {
        double temp = 0.0;
        for (int i = 0; i < numberOfAssociationSites; i++) {
            temp += dFCPAdNdXi(i, phase) * getXsitedV()[i];
        }
        return temp;
    }

    /** Getter for property xsite.
     * @return Value of property xsite.
     */
    public double[] getXsite() {
        return this.xsite;
    }

    /** Setter for property xsite.
     * @param xsite New value of property xsite.
     */
    public void setXsite(double[] xsite) {
        this.xsite = xsite;
    }

    public void setXsite(int i, double xsite) {
        this.xsite[i] = xsite;
    }

    public double[] getXsitedV() {
        return this.xsitedV;
    }

    public void setXsitedV(int i, double xsitedV) {
        this.xsitedV[i] = xsitedV;
    }

    public double[] getXsiteOld() {
        return this.xsiteOld;
    }

    /** Setter for property xsite.
     * @param xsite New value of property xsite.
     */
    public void setXsiteOld(double[] xsiteOld) {
        this.xsiteOld = xsiteOld;
    }

    public void setXsiteOld(int i, double xsiteOld) {
        this.xsiteOld[i] = xsiteOld;
    }

    public double[] getXsitedT() {
        return this.xsitedT;
    }

    public void setXsitedT(int i, double xsitedT) {
        this.xsitedT[i] = xsitedT;
    }

    public double[] getXsitedTdT() {
        return null;
    }

    public void setXsitedTdT(int i, double xsitedT) {
    }

    public void setXsitedni(int xnumb, int compnumb, double val) {
        xsitedni[xnumb][compnumb] = val;
    }
}
