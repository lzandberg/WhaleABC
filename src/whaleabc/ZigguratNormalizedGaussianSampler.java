
package whaleabc;



/*
002 * Licensed to the Apache Software Foundation (ASF) under one or more
003 * contributor license agreements.  See the NOTICE file distributed with
004 * this work for additional information regarding copyright ownership.
005 * The ASF licenses this file to You under the Apache License, Version 2.0
006 * (the "License"); you may not use this file except in compliance with
007 * the License.  You may obtain a copy of the License at
008 *
009 *      http://www.apache.org/licenses/LICENSE-2.0
010 *
011 * Unless required by applicable law or agreed to in writing, software
012 * distributed under the License is distributed on an "AS IS" BASIS,
013 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
014 * See the License for the specific language governing permissions and
015 * limitations under the License.
016 */
/**
023 * <a href="https://en.wikipedia.org/wiki/Ziggurat_algorithm">
024 * Marsaglia and Tsang "Ziggurat" method</a> for sampling from a Gaussian
025 * distribution with mean 0 and standard deviation 1.
026 *
027 * <p>The algorithm is explained in this
028 * <a href="http://www.jstatsoft.org/article/view/v005i08/ziggurat.pdf">paper</a>
029 * and this implementation has been adapted from the C code provided therein.</p>
030 *
031 * <p>Sampling uses:</p>
032 *
033 * <ul>
034 *   <li>{@link UniformRandomProvider#nextLong()}
035 *   <li>{@link UniformRandomProvider#nextDouble()}
036 * </ul>
037 *
038 * @since 1.1
039 */
public class ZigguratNormalizedGaussianSampler {
    /** Start of tail. */
    private static final double R = 3.442619855899;
    /** Inverse of R. */
    private static final double ONE_OVER_R = 1 / R;
    /** Rectangle area. */
    private static final double V = 9.91256303526217e-3;
    /** 2^63. */
    private static final double MAX = Math.pow(2, 63);
    /** 2^-63. */
    private static final double ONE_OVER_MAX = 1d / MAX;
    /** Number of entries. */
    private static final int LEN = 128;
    /** Index of last entry. */
    private static final int LAST = LEN - 1;
    /** Auxiliary table. */
    private static final long[] K = new long[LEN];
    /** Auxiliary table. */
    private static final double[] W = new double[LEN];
    /** Auxiliary table. */
    private static final double[] F = new double[LEN];
    /** Underlying source of randomness. */
    WhaleParameters rng;

    static {
        // Filling the tables.

        double d = R;
        double t = d;
        double fd = gauss(d);
        final double q = V / fd;

        K[0] = (long) ((d / q) * MAX);
        K[1] = 0;

        W[0] = q * ONE_OVER_MAX;
        W[LAST] = d * ONE_OVER_MAX;

        F[0] = 1;
        F[LAST] = fd;

        for (int i = LAST - 1; i >= 1; i--) {
            d = Math.sqrt(-2 * Math.log(V / d + fd));
            fd = gauss(d);

            K[i + 1] = (long) ((d / t) * MAX);
            t = d;

            F[i] = fd;

            W[i] = d * ONE_OVER_MAX;
        }
    }

    /**
096     * @param rng Generator of uniformly distributed random numbers.
097     */
    public ZigguratNormalizedGaussianSampler(WhaleParameters rng) {
        this.rng = rng;
    }

    /** {@inheritDoc} */

    public double sample() {
        final long j = rng.nextLong();
        final int i = (int) (j & LAST);
        //if ((j<K[i])&&(K[i]+j>0)){
        if (Math.abs(j) < K[i]) {
            return j * W[i];
        } else {
            return fix(j, i);
        }
    }


    /**
121     * Gets the value from the tail of the distribution.
122     *
123     * @param hz Start random integer.
124     * @param iz Index of cell corresponding to {@code hz}.
125     * @return the requested random value.
126     */
    private double fix(long hz,
                       int iz) {
        double x;
        double y;

        x = hz * W[iz];
        if (iz == 0) {
            // Base strip.
            // This branch is called about 5.7624515E-4 times per sample.
            do {
                y = -Math.log(rng.nextDouble());
                x = -Math.log(rng.nextDouble()) * ONE_OVER_R;
            } while (y + y < x * x);

            final double out = R + x;
            return hz > 0 ? out : -out;
        } else {
            // Wedge of other strips.
            // This branch is called about 0.027323 times per sample.
            if (F[iz] + rng.nextDouble() * (F[iz - 1] - F[iz]) < gauss(x)) {
                return x;
            } else {
                // Try again.
                // This branch is called about 0.012362 times per sample.
                return sample();
            }
        }
    }

    /**
157     * @param x Argument.
158     * @return \( e^{-\frac{x^2}{2}} \)
159     */
    private static double gauss(double x) {
        return Math.exp(-0.5 * x * x);
    }
}



























































