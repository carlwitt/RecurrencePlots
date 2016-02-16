package de.uberspace.wittcarl.phasespace;

import java.util.Arrays;

/**
 * Computes statistics of a trajectory's phase space points. Typical examples include the maximum distance between two phase space points, which can be used to estimate the recurrence threshold.
 */
public class PhaseSpaceDistribution {

    private final double approxDiameter;
    private final double[][] trajectory;
    long[] approximateDistanceDistribution;

    public PhaseSpaceDistribution(double[][] trajectory){
        this.trajectory = trajectory;
        approxDiameter = maxPhaseSpaceDiameterApproximate(this.trajectory);
    }

    public long[] approximateDistanceDistribution(){
        if(this.approximateDistanceDistribution == null) approximateDistanceDistribution = approximateDistanceDistribution(trajectory);
        return approximateDistanceDistribution;
    }

    public static double euclideanSquaredDistance(double[][] ts, int i, int j){
        double squared_distance = 0;
        int dim = ts.length;
        for (int d = 0; d < dim; d++) {
            squared_distance += Math.pow(ts[d][i] - ts[d][j], 2);// * (ts[d][col] - ts[d][row]);
        }
        return squared_distance;
    }

    public static long[] approximateSubsequentDistanceDistribution(double[][] ts, int numBins, double maxDistance){
        long distanceCount[] = new long[numBins+1];
        final double scale = numBins / maxDistance;

        int N = ts[0].length;
        for (int i = 0; i < N-1; i++) {
            double squared_distance = euclideanSquaredDistance(ts, i, i+1);
            int binIdx = (int) (scale * Math.sqrt(squared_distance));
            distanceCount[binIdx]++;
        }

        // remove trailing zeros
        int firstZeroOfTrail = distanceCount.length; // the index of the first of the trailing zeros
        for (int i = distanceCount.length-1; i > 0; i--) {
            if(distanceCount[i] > 0) break;
            firstZeroOfTrail = i;
        }

        return Arrays.copyOf(distanceCount, firstZeroOfTrail);

    }
    public static long[] approximateDistanceDistribution(double[][] ts){
        return approximateDistanceDistribution(ts, 1500, maxPhaseSpaceDiameterApproximate(ts));
    }

    /**
     * Approximates the distribution of pairwise distances between phase space points. The approximation consists in quantizing distances.
     * Relative to the approximate maximum phase space diameter {@link #maxPhaseSpaceDiameterApproximate(double[][])}
     * @param ts The trajectory [dimension][time step]
     * @param numBins The number of different distances.
     * @param upTo The maximum distance to record.
     */
    public static long[] approximateDistanceDistribution(double[][] ts, int numBins, double upTo){
        final double scale = numBins / upTo;
        long[] distanceCount = new long[numBins+1];
        final int N = ts[0].length;

        for (int col = 0; col < N; col++) {
            for (int row = 0; row < N; row++) {
                // compute euclidean distance in dim dimensions
                double squared_distance = euclideanSquaredDistance(ts, col, row);
                int binIdx = (int) (scale * Math.sqrt(squared_distance));
                if(binIdx <= numBins) distanceCount[binIdx]++;
            }
        }
//        distanceCount[0] = Double.doubleToLongBits(upTo);
        return distanceCount;
    }

    /**
     * Computes the maximum euclidean distance between any point on the first and any point on the second trajectory.
     * Quadratic runtime.
     * @param trajectory1 First trajectory. The array dimensions correspond to [phase space dimension][time step].
     * @param trajectory2 Second trajectory.
     * @return The maximum distance between any two points.
     */
    public static double maxPhaseSpaceDiameter(double[][] trajectory1, double[][] trajectory2){
        final int N = trajectory1[0].length;
        double max_squared_distance = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {
                double squared_distance = 0;
                for (int d = 0; d < trajectory1.length; d++) squared_distance += (trajectory1[d][i] - trajectory2[d][j]) * (trajectory1[d][i] - trajectory2[d][j]);
                max_squared_distance = Math.max(max_squared_distance, squared_distance);
            }
        }
        return Math.sqrt(max_squared_distance);
    }

    /**
     * Approximates the phase space diamter by computing the diagonal of the bounding (axis parallel) hypercube.
     * This is not a particularly tight bound. Better ones might be found.
     * Can be used as an auxiliary method to quickly find an upper to use e.g. for allocating a histogram array prior to counting the quantized pairwise distances.
     * @param ts The input time series [phase space dimension][time step]
     * @return An upper bound of the maximum phase space diameter.
     */
    public static double maxPhaseSpaceDiameterApproximate(double[][] ts){

        int numDims = ts.length;
        double[] mins = new double[numDims], maxs = new double[numDims];
        Arrays.fill(mins, Double.POSITIVE_INFINITY);
        Arrays.fill(maxs, Double.NEGATIVE_INFINITY);

        // find the minimum and maximum in each dimension
//        double[] mean = new double[numDims];
        for (int dim = 0; dim < numDims; dim++) {
            for (int t = 0; t < ts[dim].length; t++) {
                if(ts[dim][t] < mins[dim]) mins[dim] = ts[dim][t];
                if(ts[dim][t] > maxs[dim]) maxs[dim] = ts[dim][t];
//                mean[dim] += ts[dim][t];
            }
//            mean[dim] /= ts[dim].length;
        }


        // The length of all diagonals of the hypercube is the same (because one needs to make the step in each dimension to reach the opposite point).
        // The square root of the sum of squared side lengths.
        double sum = 0;
        for (int dim = 0; dim < numDims; dim++) {
            sum += (maxs[dim] - mins[dim]) * (maxs[dim] - mins[dim]);
        }
        return Math.sqrt(sum);

    }
}