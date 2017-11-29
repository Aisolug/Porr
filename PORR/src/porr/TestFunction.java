package porr;

import static java.lang.Math.pow;

/**
 *
 * @author aisolug
 */
public class TestFunction {
    
    int exponent;
    int[] coefficiency;
    int[][] multiplyX;
    int n;
    int initial;
    int[][] gradientX;
    int gradientExponent;
    double[] jacobiResults;
    int[][] rosenbrockX;
    double[] rosenbrockFun;
    double[][] rosenbrockInverseDM;
    double[] rosenbrockResults;
    
    
    public TestFunction(int[] coefficiency,int exponent, int x, int n){
        
        this.coefficiency = coefficiency;
        this.exponent = exponent;
        multiplyX = new int[n][n];
        this.n = n;
        initial = x;   
    }
    
    public void countFunction(){
        
        for (int j = initial-1; j< n; j++){
            for (int i = 1; i < initial; i++){
                multiplyX[j-i][j-i] += coefficiency[initial-1-i];
           }
            multiplyX[j][j] += j*coefficiency[2]; 
            for (int i = 0; i < initial-1; i++){
                multiplyX[j-2][j-i] += (j+1)*coefficiency[initial+initial-2-i];
            }  
        }
        multiplyX[0][0] += 5; //zapewniona dominacja
    }
    
    public void complementMatrix(){
        
        for (int i = 1; i < n; i ++){
            for (int k = 0; k < n; k ++){
                if(k < i){
                    multiplyX[i][k] = multiplyX[k][i];    
                }
            }
        }
    }
    
    public void derivativeMatrix(){
        
        gradientX = multiplyX;
        for (int i = 0; i < n; i++){
            gradientX[i][i] = multiplyX[i][i]*exponent;
        }
        gradientExponent = exponent-1; 
    }
    
    public void rosenbrockDMatrix(int x0){
        rosenbrockX = new int[n][n];
        for (int i = 0; i < n; i++){
            rosenbrockX[i][i] = gradientX[i][i]*gradientExponent*(int)pow(x0,(gradientExponent-1));
        }
    }
    
    public void inverseRosenbrockDMatrix(){
        rosenbrockInverseDM = new double[n][n];
        for (int i = 0; i < n; i ++){
            rosenbrockInverseDM[i][i] = 1/((double) rosenbrockX[i][i]);
        }
    }
    
    public void rosenbrockFunction(int x0){
        
        rosenbrockFun = new double[n];
        for (int i =0; i < n; i ++){
            for (int k = 0; k < n; k++){
                if(k != i) rosenbrockFun[i] += gradientX[i][k]*x0;
                else rosenbrockFun[i] += gradientX[i][i]*pow(x0, gradientExponent);
            }
        }
    }
    
    public void resultJacobi(int x0){
        
        jacobiResults = new double[n];
        double[] sum = new double[n];
        double divider = (double) gradientExponent;
        for (int i = 0; i < n; i ++){
            for (int k = 0; k < n; k ++){
                if(i != k){
                    sum[i] -= gradientX[i][k]*x0;
                }
            }
        }
        for (int i = 0; i < n; i++){
            
            jacobiResults[i] = pow(Math.E, Math.log(sum[i]/gradientX[i][i])/divider);
        }
    }
    
    public double[] getJacobiResults(){
        
        return jacobiResults;
    }

    public void RosenbrockEquation(int x0){
        
        rosenbrockResults = new double[n];
        for (int i = 0; i < n; i ++){
            rosenbrockResults[i] = x0 - rosenbrockInverseDM[i][i]*rosenbrockFun[i];
        }  
    }
    
    public void printFunction(){
        
        for (int i = 0; i < n; i ++){
            for (int k = 0; k < n; k ++){
                if (multiplyX[i][k] != 0){    
                    if (k == i) System.out.print("+ "+multiplyX[i][i]+"x["+(i+1)+"]^"+exponent);
                    else System.out.print(" "+multiplyX[i][k]+"x["+(i+1)+"]*"+"x["+(k+1)+"]");
                }
            }
        }
        System.out.println("\n");
    }
    
    public void printJacobianEquations(){
        
        for (int k = 0; k < n; k ++){
            for (int i = 0; i < n; i ++){
                if (gradientX[i][k] != 0){    
                    if (k == i) System.out.print("+ "+multiplyX[i][i]+"x["+(i+1)+"]^"+gradientExponent);
                    else System.out.print(" "+multiplyX[i][k]+"x["+(i+1)+"]");
                }  
            }
            System.out.print(" = 0\n");
        }  
        System.out.println();
    }
    
    public void printJacobiResults(){
        for(int i = 0; i < n; i++){
            System.out.println(jacobiResults[i]);
        }
    }
    
    public void printRosenbrockFunction(){
        for (int i = 0; i < n; i++){
            System.out.println(rosenbrockFun[i]);
        }
    }
    
    public void printRosenbrockEquation(){
        
        for (int i = 0; i < n; i ++){
            System.out.println(rosenbrockResults[i]);
        }  
    }
    
    public void printMatrix(){
        
        for (int i = 0; i < n; i ++){
            for (int k = 0; k < n; k ++){
                System.out.print(multiplyX[i][k]+" ");
            }
            System.out.println();   
        }
        System.out.println();
    }
    
    public void printDerivativeMatrix(){
        
        for (int i = 0; i < n; i ++){
            for (int k = 0; k < n; k ++){
                System.out.print(gradientX[i][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }
        
    public void printRosenbrockMatrix(){
        
        for (int i = 0; i < n; i ++){
            for (int k = 0; k < n; k ++){
                System.out.print(rosenbrockX[i][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }
    
     public void printInverseRosenbrockDM(){
         
         for (int i = 0; i < n; i ++){
            for (int k = 0; k < n; k ++){
                System.out.print(rosenbrockInverseDM[i][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
     }
}