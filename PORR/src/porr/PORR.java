package porr;

/**
 *
 * @author aisolug
 */
public class PORR {


    //int[] function = new int[2];
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        int[] coefficiency = {7,20,5,-9,-10};
        int exponent = 4;
        int x0 = 2;
        int n = 50;

        TestFunction test = new TestFunction(coefficiency, exponent, 3, n);
        
        test.countFunction();        
        test.printFunction();
        
        test.complementMatrix();
        test.printMatrix();
        
        test.derivativeMatrix();
        test.printDerivativeMatrix();
        test.printJacobianEquations();
        test.resultJacobi(x0); //podstawienie za x0
        test.printJacobiResults();
        System.out.println();
        
        test.rosenbrockDMatrix(x0); //podstawienie za x0
        //test.printRosenbrockMatrix();
        
        test.rosenbrockFunction(x0); //podstawienie za x0
        //test.printRosenbrockFunction();
        
        test.inverseRosenbrockDMatrix();
        //test.printInverseRosenbrockDM();
        
        test.RosenbrockEquation(x0);
        test.printRosenbrockEquation();
        
        
        
        
    }
    
}
