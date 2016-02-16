package de.uberspace.wittcarl.experimental;

import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import de.uberspace.wittcarl.phasespace.PhaseSpaceDistribution;
import junit.framework.TestCase;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 20.08.15.
 */
public class DimensinoalityTest extends TestCase {

    public void testPhaseSpaceTrajectoryDiameterExpansionRates() throws IOException {

        // generate a base trajectory of the standard lorenz system
        final int numSteps = 5000;

        String outDir = String.format("data/junit-test-data/noisy-%s/", numSteps);
//        Files.createDirectories(Paths.get(outDir));

        TimeSeriesGenerator timeSeriesGenerator = new TimeSeriesGenerator(30, numSteps, null);

        double[][][][] allTrajectories = timeSeriesGenerator.getAllTrajectories();

        for (int i = 0; i < allTrajectories.length; i++) {
            String name = TimeSeriesGenerator.system_names[i];
            System.out.println(name);
            DescriptiveStatistics statistics = new DescriptiveStatistics();

            for (int j = 0; j < allTrajectories[i].length; j++) {
                double[][] original = allTrajectories[i][j];

                ArrayList<Double> diameters = new ArrayList<Double>();
                // add noise and write the trajectories to disk


            /*for (double noiseRatio : DRQA_All.noiseRatiosById.values()) {*/
                for (double noiseRatio : new double[]{0.0, 0.8}) {
                    double[][] noisyTrajectory = TimeSeriesGenerator.addNoise(original, 0.05, noiseRatio);
                    String filename = String.format("pst_%s-%s-%s.csv", name, numSteps, noiseRatio);
//                DRQA drqa = new DRQA(noisyTrajectory, noisyTrajectory, 0.05);
                    diameters.add(PhaseSpaceDistribution.maxPhaseSpaceDiameter(noisyTrajectory, noisyTrajectory));
//                Files.write(Paths.get(outDir + filename), TimeSeriesGenerator.trajectoryToCsvString(noisyTrajectory).toString().getBytes());
                }

                final Double first = diameters.get(0);
                final Double last = diameters.get(diameters.size() - 1);
                statistics.addValue(last/first);
//                System.out.println(first);
//                System.out.println(last);
//                System.out.println((last / first));
//                System.out.println("\n");
            }
            System.out.println(String.format("diamter expansion Min: %s", statistics.getMin()));
            System.out.println(String.format("diamter expansion Mean: %s", statistics.getMean()));
            System.out.println(String.format("diamter expansion Max: %s\n", statistics.getMax()));


        }
    }

}
