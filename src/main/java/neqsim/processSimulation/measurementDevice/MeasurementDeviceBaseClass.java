/*
 * ProcessEquipmentBaseClass.java
 *
 * Created on 6. juni 2006, 15:12
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package neqsim.processSimulation.measurementDevice;

import neqsim.processSimulation.measurementDevice.online.OnlineSignal;

/**
 *
 * @author ESOL
 */
public abstract class MeasurementDeviceBaseClass implements MeasurementDeviceInterface {

    private static final long serialVersionUID = 1000;

    /**
     * 
     * @return the onlineSignal
     */
    public OnlineSignal getOnlineSignal() {
        return onlineSignal;
    }

    /**
     * @param onlineSignal the onlineSignal to set
     */
    public void setOnlineSignal(OnlineSignal onlineSignal) {
        this.onlineSignal = onlineSignal;
    }

    /**
     * @return the isOnlineSignal
     */
    public boolean isOnlineSignal() {
        return isOnlineSignal;
    }

    /**
     * @param isOnlineSignal the isOnlineSignal to set
     */
    public void setIsOnlineSignal(boolean isOnlineSignal, String plantName, String transmitterame) {
        this.isOnlineSignal = isOnlineSignal;
        onlineSignal = new OnlineSignal(plantName, transmitterame);
    }

    protected String name = "default";
    protected String unit = "-";
    private double maximumValue = 1.0;
    private double minimumValue = 0.0;
    private boolean logging = false;
    private OnlineSignal onlineSignal = null;
    private boolean isOnlineSignal = false;

    /**
     * Creates a new instance of ProcessEquipmentBaseClass
     */
    public MeasurementDeviceBaseClass() {
    }

    public void displayResult() {
    }

    ;
    public String getName() {
        return name;
    }

    public String getUnit() {
        return unit;
    }

    public void setName(String nameset) {
        name = nameset;
    }

    public double getMeasuredValue() {
        return 0.0;
    }

    public void setUnit(String unit) {
        this.unit = unit;
    }

    public double getMaximumValue() {
        return maximumValue;
    }

    public void setMaximumValue(double maximumValue) {
        this.maximumValue = maximumValue;
    }

    public double getMinimumValue() {
        return minimumValue;
    }

    public void setMinimumValue(double minimumValue) {
        this.minimumValue = minimumValue;
    }

    public double getMeasuredPercentValue() {
        return (getMeasuredValue() - minimumValue) / (maximumValue - minimumValue) * 100;
    }

    public boolean isLogging() {
        return logging;
    }

    public void setLogging(boolean logging) {
        this.logging = logging;
    }

    public double getOnlineValue() {
        return getOnlineSignal().getValue();
    }
}
