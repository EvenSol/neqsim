/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package neqsim.processSimulation.processEquipment.util;

import neqsim.processSimulation.processEquipment.ProcessEquipmentBaseClass;
import neqsim.processSimulation.processEquipment.stream.Stream;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicOperations.ThermodynamicOperations;

/**
 *
 * @author esol
 */
public class StreamSaturatorUtil extends ProcessEquipmentBaseClass {

    private static final long serialVersionUID = 1000;

    Stream inletStream;
    Stream outStream;
    SystemInterface thermoSystem;
    ThermodynamicOperations thermoOps;

    public StreamSaturatorUtil(Stream inletStream) {
        setInletStream(inletStream);
    }

    public void setInletStream(Stream inletStream) {
        this.inletStream = inletStream;

        thermoSystem = (SystemInterface) inletStream.getThermoSystem().clone();
        outStream = new Stream(thermoSystem);
    }

    public Stream getOutStream() {
        return outStream;
    }

    public void run() {
        System.out.println("valve running..");
        thermoSystem = (SystemInterface) inletStream.getThermoSystem().clone();
        thermoOps = new ThermodynamicOperations(thermoSystem);
        thermoSystem.init(3);

        outStream.setThermoSystem(thermoSystem);
    }
}
