import java.util.Arrays;

/**
 * Created by Jordan Crawford on 9/13/17.
 *
 *
 * Sample output:

 Problem 1: Original Matrix:
 [1.00,3.00,2.00,1.00,-2.00,]
 [4.00,2.00,1.00,2.00,2.00,]
 [2.00,1.00,2.00,3.00,1.00,]
 [1.00,2.00,4.00,1.00,-1.00,]
 After Elimination:
 [1.00,3.00,2.00,1.00,-2.00,]
 [0.00,-10.00,-7.00,-2.00,10.00,]
 [0.00,0.00,1.50,2.00,0.00,]
 [0.00,0.00,0.00,-3.40,0.00,]
 Result of backsub:
 [1.0, -1.0, 0.0, -0.0]


 ===Problem 2===
 *Delta: 0.06250000000000000000
 Actual solution:     x1 = 1.0666666667 x2 = 0.9333333333
 Calculated solution: x1 = 1.0666666667 x2 = 0.9333333333
 *Delta: 0.00390625000000000000
 Actual solution:     x1 = 1.0039215686 x2 = 0.9960784314
 Calculated solution: x1 = 1.0039215686 x2 = 0.9960784314
 *Delta: 0.00024414062500000000
 Actual solution:     x1 = 1.0002442002 x2 = 0.9997557998
 Calculated solution: x1 = 1.0002442002 x2 = 0.9997557998
 *Delta: 0.00001525878906250000
 Actual solution:     x1 = 1.0000152590 x2 = 0.9999847410
 Calculated solution: x1 = 1.0000152590 x2 = 0.9999847410
 *Delta: 0.00000095367431640625
 Actual solution:     x1 = 1.0000009537 x2 = 0.9999990463
 Calculated solution: x1 = 1.0000009537 x2 = 0.9999990463
 *Delta: 0.00000005960464477539
 Actual solution:     x1 = 1.0000000596 x2 = 0.9999999404
 Calculated solution: x1 = 1.0000000596 x2 = 0.9999999404
 *Delta: 0.00000000372529029846
 Actual solution:     x1 = 1.0000000037 x2 = 0.9999999963
 Calculated solution: x1 = 1.0000000000 x2 = 0.9999999963
 *Delta: 0.00000000023283064365
 Actual solution:     x1 = 1.0000000002 x2 = 0.9999999998
 Calculated solution: x1 = 1.0000000000 x2 = 0.9999999998
 *Delta: 0.00000000001455191523
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 1.0000000000 x2 = 1.0000000000
 *Delta: 0.00000000000090949470
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 1.0000000000 x2 = 1.0000000000
 *Delta: 0.00000000000005684342
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 1.0000000000 x2 = 1.0000000000
 *Delta: 0.00000000000000355271
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 1.0000000000 x2 = 1.0000000000
 *Delta: 0.00000000000000022204
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 1.0000000000 x2 = 1.0000000000
 *Delta: 0.00000000000000001388
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 0.0000000000 x2 = 1.0000000000
 *Delta: 0.00000000000000000087
 Actual solution:     x1 = 1.0000000000 x2 = 1.0000000000
 Calculated solution: x1 = 0.0000000000 x2 = 1.0000000000

 */
public class Matrix {

    double[][] data;

    public Matrix( double[][] data ) {
        this.data = data;
    }

    public void swap(int row1, int row2) {
        double[] temp = data[row1];
        data[row1] = data[row2];
        data[row2] = temp;
    }


    public void naiveGaussianEliminate() {
        for(int row = 0; row < data.length; row ++) {
            double pivot = data[row][row];
            if( pivot == 0 ) {
                for( int row_prime = row + 1; row_prime < data.length; row_prime ++ ) {
                    if( data[row_prime][row] != 0 ) {
                        swap(row, row_prime);
                        break;
                    }
                }
                if( pivot == 0 ) {
                    throw new RuntimeException("Nonsingular matrix!");
                }
            }
            for( int row_prime = row + 1; row_prime < data.length; row_prime ++) {
                double factor = data[row_prime][row]/pivot;
                for( int col = row; col < data[row_prime].length; col++ ) {
                    data[row_prime][col] = data[row_prime][col] - factor * data[row][col];
                }
            }
        }
    }

    public double[] backsub() {
        double[] output = new double[data.length];
        for( int row = data.length - 1; row >=0; row --) {
            double result = data[row][data[row].length - 1];
            for( int col = data.length - 1; col >= row; col --) {
                result -= data[row][col] * output[col];
            }
            output[row] = result/data[row][row];
        }
        return output;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for(int i = 0; i < data.length; i ++ ) {
            for(int j = 0; j < data[i].length; j++) {
                sb.append(String.format("%3.2f,", data[i][j] ) );
            }
            sb.append("]");
            if( i +1 != data.length) {
                sb.append("\n[");
            }
        }
        return sb.toString();
    }

    public static void main( String[] args ) {
        double[][] data = new double[][] {
                new double[] {1,3,2,1,-2},
                new double[] {4,2,1,2,2},
                new double[] {2,1,2,3,1},
                new double[] {1,2,4,1,-1}};

        Matrix mat = new Matrix( data );

        System.out.printf("Problem 1: Original Matrix:\n%s\n", mat);

        mat.naiveGaussianEliminate();

        System.out.printf("After Elimination:\n%s\n", mat);

        double[] backsubResult = mat.backsub();

        System.out.printf("Result of backsub:\n%s\n", Arrays.toString( backsubResult ) );

        System.out.println("\n\n===Problem 2===");
        for(int i = 1; i < 16; i ++) {
            double delta = Math.pow(2, -4*i);
            System.out.printf("*Delta: %.20f\n", delta);
            System.out.printf("Actual solution:     x1 = %2.10f x2 = %2.10f\n", 1/(1-delta), (1-2* delta)/(1-delta) );
            Matrix toSolve = new Matrix(new double[][] {
                    new double[] {delta, 1, 1},
                    new double[] {1, 1, 2}
            });
            toSolve.naiveGaussianEliminate();
            double[] result = toSolve.backsub();
            System.out.printf("Calculated solution: x1 = %2.10f x2 = %2.10f\n", result[0], result[1] );

        }
    }
}
