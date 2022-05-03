package com.LyapXool;

import java.util.Arrays;

public class ArrayOperations {
    static public double dot(double[] vA, double[] vB)
    {
        int sa = vA.length;
        int sb = vB.length;
        if (sa != sb)
        {
            System.out.println("The vectors must have the same size!");
            System.exit(0);
        }
        double product = 0.0;

        // Loop for calculate cot product
        for (int i = 0; i < sa; i++)
            product += vA[i] * vB[i];
        return product;
    }
    static public double[] vdiff(double[] vA, double[] vB)
    {
        int sa = vA.length;
        int sb = vB.length;
        if (sa != sb)
        {
            System.out.println("The vectors must have the same size!");
            System.exit(0);
        }
        double[] diff = new double[sa];
        Arrays.fill(diff, 0.0);

        for (int i = 0; i < sa; i++)
            diff[i] = vA[i] - vB[i];
        return diff;
    }
    static public void FillMatrixConst(double[][] matrixToFill, double constVal)
    {
        int minRows = matrixToFill.length;
        int minCols = matrixToFill[0].length;
        for (int i = 0; i < minRows; ++i)
        {
            for (int j = 0; j < minCols; ++j)
            {
                matrixToFill[i][j] = constVal;
            }
        }
    }

    static double[][] ResizeArray(double[][] matrix, int rows, int cols) {
        double[][] temp = new double[rows][cols];
        for (int i = 0; i < rows; i++)
            for (int j=0; j<cols;++j)
                temp[i][j]=0.0;
        return temp;
    }


    private static final double EPSILON = 1e-10;

    // is symmetric
    public static boolean isSymmetric(double[][] A) {
        int N = A.length;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                if (Math.abs(A[i][j] - A[j][i])>10e-10) return false;
            }
        }
        return true;
    }

    // is symmetric
    public static boolean isSquare(double[][] A) {
        int N = A.length;
        for (int i = 0; i < N; i++) {
            if (A[i].length != N) return false;
        }
        return true;
    }
    public static double[][] transpose(double[][] matrixToTranspose){
        int rows=matrixToTranspose.length;
        int cols=matrixToTranspose[0].length;
        double[][] transposed=new double[rows][cols];
        for(int i=0; i<rows; ++i)
        {
            for(int j=0; j<cols; ++j)
            {
                transposed[j][i]=matrixToTranspose[i][j];
            }
        }
        return transposed;
    }
    public static double[] multiply(double[][] a, double[] x) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += a[i][j] * x[j];
        return y;
    }
    public static double[][] invert(double a[][])
    {
        int n = a.length;
        double x[][] = new double[n][n];
        double b[][] = new double[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i)
            b[i][i] = 1;

        // Transform the matrix into an upper triangle
        gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                            -= a[index[j]][i]*b[index[i]][k];

        // Perform backward substitutions
        for (int i=0; i<n; ++i)
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j)
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k)
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }
    public static void gaussian(double a[][], int index[])
    {
        int n = index.length;
        double c[] = new double[n];

        // Initialize the index
        for (int i=0; i<n; ++i)
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i)
        {
            double c1 = 0;
            for (int j=0; j<n; ++j)
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j)
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i)
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1)
                {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i)
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }
    public static double[] vectorMultiply(double s, double[] v) {
        double[] result = new double[v.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = s * v[i];
        }
        return result;
    }
}

