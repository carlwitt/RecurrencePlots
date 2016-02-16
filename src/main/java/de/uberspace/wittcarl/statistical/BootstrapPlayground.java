package de.uberspace.wittcarl.statistical;

import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 23.10.15.
 *
 * Experiments with bootstrapping line length histograms (to obtain confidence intervals for RQA measures).
 *
 */
public class BootstrapPlayground {

    public void testHistogramBootstrapping(){

        final int resamples = 100000;
        final double p_value = 0.01;

        TimeSeriesGenerator tsg = new TimeSeriesGenerator(1, 5000, "lorenz");
        double[][] trajectory = tsg.lorenz_standard[0];
        DRQA drqa = new DRQA(trajectory, trajectory, 0.05);
        drqa.computeRQA(2, 2, 2);

        double[] means = new double[4], cis = new double[4];
        System.out.println(String.format("Bootstrapping sampling distribution, %s resamples, p-value=%s", resamples, p_value));
        DRQA.bootstrapSamplingDistribution(p_value, resamples, 2, drqa.l_hist, means, cis);

        String[] names = new String[]{"Filter Ratio", "Mean", "Median", "Entropy"};
        DRQA.HistogramStatistic[] statistics = new DRQA.HistogramStatistic[]{DRQA.HistogramStatistic.FILTER_RATIO, DRQA.HistogramStatistic.AVERAGE, DRQA.HistogramStatistic.MEDIAN, DRQA.HistogramStatistic.ENTROPY};
        for (int i = 0; i < 4; i++) {
            System.out.println(names[i]);
            System.out.println(String.format("Original Value: %s", drqa.diagonal_rqa_definite.get(statistics[i])));
            System.out.println(String.format("Bootstrap Mean: %s", means[i]));
            System.out.println(String.format("Bootstrap :CI width: %s", cis[i]));
            System.out.println("");
        }
    }

}
