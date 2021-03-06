/*
 * Heater.java
 *
 * Created on 15. mars 2001, 14:17
 */

package neqsim.processSimulation.processEquipment.heatExchanger;

import neqsim.processSimulation.processEquipment.ProcessEquipmentBaseClass;
import neqsim.processSimulation.processEquipment.ProcessEquipmentInterface;
import neqsim.processSimulation.processEquipment.stream.StreamInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicOperations.ThermodynamicOperations;

/**
 *
 * @author  Even Solbraa
 * @version
 */
public class ReBoiler extends ProcessEquipmentBaseClass implements ProcessEquipmentInterface{

    private static final long serialVersionUID = 1000;
    
    ThermodynamicOperations testOps;
    boolean setTemperature=false;
    String name= new String();
    StreamInterface outStream;
    StreamInterface inStream;
    SystemInterface system;
    private double reboilerDuty = 0.0;
    /** Creates new Heater */
    public ReBoiler() {
    }
    
    public ReBoiler(StreamInterface inStream) {
        this.inStream = inStream;
        outStream = (StreamInterface) inStream.clone();
    }
    
    public void setName(String name){
        this.name = name;
    }
    
    
    public StreamInterface getOutStream(){
        return outStream;
    }
    
    public void run(){
        system = (SystemInterface) inStream.getThermoSystem().clone();
        testOps = new ThermodynamicOperations(system);
        testOps.TPflash();
        double oldH = system.getEnthalpy();
        testOps = new ThermodynamicOperations(system);
        testOps.TPflash();
        testOps.PHflash(oldH+reboilerDuty,0);
        outStream.setThermoSystem(system);
//        if(setTemperature) system.setTemperature(temperatureOut);
//        else system.setTemperature(system.getTemperature()+dT);
//        testOps = new ThermodynamicOperations(system);
//        system.setTemperat ure(temperatureOut);
//        testOps.TPflash();
//        double newH = system.getEnthalpy();
//        dH = newH - oldH;
//        // system.setTemperature(temperatureOut);
//        //  testOps.TPflash();
//        //    system.setTemperature(temperatureOut);
//        outStream.setThermoSystem(system);
    }
    
    public void displayResult(){
       System.out.println("out Temperature " + reboilerDuty);
    }
    
    public String getName() {
        return name;
    }
    
    public void runTransient() {
    }
    
    public double getReboilerDuty() {
        return reboilerDuty;
    }
    
    public void setReboilerDuty(double reboilerDuty) {
        this.reboilerDuty = reboilerDuty;
    }
    
}
