/*
 * LevenbergMarquardt.java
 *
 * Created on 22. januar 2001, 23:00
 */

package neqsim.statistics.parameterFitting.nonLinearParameterFitting;


import Jama.*;
import static cern.jet.stat.Gamma.beta;
import static cern.jet.stat.Probability.chiSquare;
import neqsim.statistics.parameterFitting.StatisticsBaseClass;
/**
 *
 * @author  Even Solbraa
 * @version
 */
public class LevenbergMarquardt extends StatisticsBaseClass {

    private static final long serialVersionUID = 1000;
    double oldChiSquare=1e100;
    double newChiSquare=0;
    Matrix parameterStdDevMatrix;
    Matrix parameterUncertaintyMatrix;
    Thread thisThread;
    boolean solved=false;
    private int maxNumberOfIterations = 50;
    
    
    /** Creates new LevenbergMarquardt */
    public LevenbergMarquardt() {
        thisThread = new Thread();
    }
    
    public Object clone(){
        LevenbergMarquardt clonedClass = null;
        try{
            clonedClass = ( LevenbergMarquardt) super.clone();
        }
        catch(Exception e) {
            e.printStackTrace(System.err);
        }
        
        return clonedClass;
    }
    
    
    public void init(){
        chiSquare = calcChiSquare();
        System.out.println("Chi square: "+ chiSquare);
        dyda = calcDerivatives();
        beta = calcBetaMatrix();
        alpha = calcAlphaMatrix();
    }
    
    public void solve(){
        setFittingParameters(sampleSet.getSample(0).getFunction().getFittingParams());
        Matrix betaMatrix;
        Matrix newParameters;
        int n=0;
        init();
        do{
            n++;
            if(chiSquare<oldChiSquare) {
                oldChiSquare = chiSquare;
            }
            betaMatrix = new Matrix(beta,1).transpose();
            Matrix alphaMatrix = new Matrix(alpha);
            betaMatrix.print(10,3);
            alphaMatrix.print(10,3);
            Matrix solvedMatrix = alphaMatrix.solve(betaMatrix);
            
            Matrix oldParameters = new Matrix(sampleSet.getSample(0).getFunction().getFittingParams(),1).copy();
            newParameters = oldParameters.copy().plus(solvedMatrix.transpose());
            Matrix diffMat = newParameters.copy().minus(oldParameters);
            this.checkBounds(newParameters);
            this.setFittingParameters(newParameters.copy().getArray()[0]);
            newChiSquare = calcChiSquare();
            if(newChiSquare>=oldChiSquare || Double.isNaN(newChiSquare)){
                newChiSquare = oldChiSquare;
                multiFactor *= 10.0;
                System.out.println("keeping old values");
                this.setFittingParameters(oldParameters.getArray()[0]);
            }
            else{
                multiFactor /= 10.0;
                System.out.println("Chi square ok : ");
            }
            System.out.println("Chi after : "+ newChiSquare);
            System.out.println("iterations: " +n);
            calcAbsDev();
            try{
                Thread.sleep(50);
            }
            catch(Exception e){}
            System.out.println("Parameters:");
            newParameters.print(100,100);
            init();
        }
        while(((Math.abs(newChiSquare)>1e-6) && n<maxNumberOfIterations && Math.abs(betaMatrix.norm2())>1.0e-6) || n<5);
        solved = true;
        
        System.out.println("iterations: " + n);
        System.out.println("Chi square: "+ chiSquare);
        System.out.println("Parameters:");
        newParameters.print(10,10);
    }
    
    
    public static void main(String[] args){
      /*  LevenbergMarquardt optim = new LevenbergMarquardt();
        TestFunction testFunction = new TestFunction();
      //  optim.setFunction(testFunction);
       
        SampleValue[] sample = new SampleValue[3];
        double sample1[] = { 6 };
        sample[0] = new SampleValue(8.5,0.1,sample1);
        double sample2[] = { 4 };
        sample[1] = new SampleValue(5.5,0.1,sample2);
        double sample3[] = { 4 };
        sample[2] = new SampleValue(5.51,0.1,sample3);
       
        SampleSet sampleSet = new SampleSet(sample);
        sampleSet = sampleSet.createNewNormalDistributedSet();
        optim.setSampleSet(sampleSet);
       // optim.solve();
        optim.runMonteCarloSimulation();*/
    }

    /**
     * @return the maxNumberOfIterations
     */
    public int getMaxNumberOfIterations() {
        return maxNumberOfIterations;
    }

    /**
     * @param maxNumberOfIterations the maxNumberOfIterations to set
     */
    public void setMaxNumberOfIterations(int maxNumberOfIterations) {
        this.maxNumberOfIterations = maxNumberOfIterations;
    }
    
}
