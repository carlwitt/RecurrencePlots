package de.uberspace.wittcarl.executable;

import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

import java.io.IOException;
import java.nio.file.*;
import java.util.Arrays;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 11.08.15.
 *
 * Used to estimate the minimum time series length for (D)RQA measures to converge.
 *
 */
public class DRQAConvergence {

    public static void main(String[] args) throws IOException {

        if(args.length < 3){
            System.err.println("Usage: max_length start_length drqa_output [p_value] [bootstrap_resamples]");
            System.exit(0);
        }

        int maxLength = Integer.parseInt(args[0]);
        int startLength = Integer.parseInt(args[1]);
        String drqaOutputFile = args[2];
        double p_value = args.length > 3 ? Double.parseDouble(args[3]) : 0.01;
        int bootstrap_resamples = args.length > 4 ? Integer.parseInt(args[4]) : 10000;

        try {
            DRQAConvergence.computeConvergenceTimes(startLength, maxLength, drqaOutputFile, p_value, bootstrap_resamples);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Generates a lorenz system time series of maxLength
     * @param initialLength the window size to start with
     * @param maxLength the maximum window size to consider (might not be included exactly)
     * @param drqaOutputFile the file path to store the DRQA measures of each iteration
     */
    public static void computeConvergenceTimes(int initialLength, int maxLength, String drqaOutputFile, double p_value, int resamples) throws IOException {

        // Generate time series of given length
        TimeSeriesGenerator timeSeriesGenerator = new TimeSeriesGenerator(1, maxLength, null);
        double[][][][] trajectories = timeSeriesGenerator.getAllTrajectories();

        final int num_measures = DRQA.LineType.values().length * DRQA.HistogramStatistic.values().length;

        // write header for DRQA output
        final Path outpath = Paths.get(drqaOutputFile);
        StringBuilder ci_header_builder = new StringBuilder();
        for (String linetype : new String[]{"L", "V", "W"})
            for (String mean_or_ci : new String[]{"MEAN","CI"})
                for (String statistic : new String[]{"FILTER_RATIO", "AVERAGE", "MEDIAN", "ENTROPY"})
                    ci_header_builder.append(linetype).append("_").append(mean_or_ci).append("_").append(statistic).append(DRQA.separator);
        ci_header_builder.deleteCharAt(ci_header_builder.length()-1);
        String header = String.format("window_size\tdrqa_computation_time\tresampling_time\t%s\t%s\n",DRQA.separatedHeader(), ci_header_builder.toString());
        Files.write(outpath, header.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);

        long started;

        for (int i = 0; i < TimeSeriesGenerator.system_names.length; i++) {

            System.out.println(TimeSeriesGenerator.system_names[i]);
            double[][] trajectory = trajectories[i][0];
            int dimensions = trajectory.length;
            double[][][] sampling = new double[][][]{new double[][]{new double[4], new double[4]}, new double[][]{new double[4], new double[4]}, new double[][]{new double[4], new double[4]}};
            int L=0, V=1, W=2, MEAN=0, CI=1;
            // increase time series length by 10% in each iteration
            for (int length = initialLength; length < maxLength; length = (int) Math.ceil(1.1*length)) {

                // extract window of according length starting at the beginning of the time series
                double[][] window = new double[dimensions][length];
                for (int dim = 0; dim < dimensions; dim++) window[dim] = Arrays.copyOfRange(trajectory[dim], 0, length);

                // compute DRQA measures
                started = System.currentTimeMillis();
                DRQA drqa = new DRQA(window, window, 0.05);
                drqa.computeRQA(2, 2, 2);

                // report time needed for current length
                final long computation_time = (System.currentTimeMillis() - started);
                System.out.println(String.format("Computed DRQA for length %s in %.2f seconds.", length, computation_time/1000.));

                started = System.currentTimeMillis();
                // compute confidence intervals
                DRQA.bootstrapSamplingDistribution(p_value, resamples, 2, drqa.l_hist, sampling[L][MEAN], sampling[L][CI]);
                DRQA.bootstrapSamplingDistribution(p_value, resamples, 2, drqa.v_hist, sampling[V][MEAN], sampling[V][CI]);
                DRQA.bootstrapSamplingDistribution(p_value, resamples, 2, drqa.w_hist, sampling[W][MEAN], sampling[W][CI]);

                // report time needed for current length
                final long resampling_time = (System.currentTimeMillis() - started);
                System.out.println(String.format("Computed Bootstrap Distributions in %s ms.", resampling_time));

                // package all results into a single string
                StringBuilder builder = new StringBuilder();
                for (int line_type = 0; line_type < 3; line_type++)
                    for (int mean_or_ci = 0; mean_or_ci < 2; mean_or_ci++)
                        for (int statistic = 0; statistic < 4; statistic++)
                            builder.append(sampling[line_type][mean_or_ci][statistic]).append(DRQA.separator);
                // remove last separator
                builder.deleteCharAt(builder.length() - 1);



                // write DRQA results
                String drqa_string = drqa.toSeparatedString("IN_MEMORY", TimeSeriesGenerator.system_names[i], DRQA.EmbeddingMethod.ORIGINAL_TRAJECTORY.ordinal(), dimensions, -1, DRQA.LMinMethod.L_MIN_FIX.ordinal(), 2, 0, length);
                String line = String.format("%s\t%s\t%s\t%s\t%s\n", length, computation_time, resampling_time, drqa_string, builder.toString());
                Files.write(outpath, line.getBytes(), StandardOpenOption.APPEND);

            } // end different window lengths
        }




    }
}
