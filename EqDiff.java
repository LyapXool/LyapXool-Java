package com.LyapXool;

public class EqDiff {
    static public void eqdff(boolean normal, double[] x, double[] f)
    {
        f[0] = -1.0 * x[0] * ((x[0] * x[0]) + (x[1] * x[1]) - (1.0 / 4.0)) * ((x[0] * x[0]) + (x[1] * x[1]) - 1.0) - x[1];
        f[1] = -1.0 * x[1] * ((x[0] * x[0]) + (x[1] * x[1]) - (1.0 / 4.0)) * ((x[0] * x[0]) + (x[1] * x[1]) - 1.0) + x[0];
        int end = x.length;
        if(normal)
        {
            double norma = Math.sqrt(ArrayOperations.dot(f,f));
            for (int i = 0; i < end; ++i) {
                f[i] /= norma;
            }

        }
    }
}
