package de.uberspace.wittcarl;

import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import de.uberspace.wittcarl.phasespace.PhaseSpaceDistribution;
import junit.framework.TestCase;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

public class TimeSeriesGeneratorTest extends TestCase{

    public void testLorenzSystem(){

        TimeSeriesGenerator tsg = new TimeSeriesGenerator(2000);

        long before = System.currentTimeMillis();
        DRQA drqa = new DRQA(tsg.lorenz_standard[0], 5);
        long after = System.currentTimeMillis();
        System.out.println(String.format("DRQA computation time: %s s", 1. * (after - before) / 1000));

    }

    public void testToString(){

        String str = DRQA.histogramToString(new long[]{0,0,7,0,0,0});
        System.out.println(str);

        int places_class = (int) Math.ceil(Math.log10(9));
        int places_timeseries = (int) Math.ceil(Math.log10(500));
        String filename = String.format(String.format("pst_%%0%sd_%%s%%0%dd.txt",places_class,places_timeseries), 0, "oscillator", 1);
        System.out.println(String.format("filename: %s", filename));

    }

    /**
     * Under very high noise levels, the RQA measures depend only on the dimensionality of the trajectory.
     * This is independent of whether the displacement vector magnitudes are corrected for dimensionality.
     * To understand this, consider for instance a 3D random walk. It has a lower chance of recurrence than a 1D random walk.
     * If the thresholds were chosen relative to the intial trajectories diameter, all should have much lower recurrence rates, but since it is chosen according to the noisy trajectories diameter, the RRs depend in this way.
     */
    public void testDimensionalityFlaw(){

        System.out.println("TimeSeriesGeneratorTest.testDimensionalityFlaw");
        long tic = System.currentTimeMillis();
        TimeSeriesGenerator tsg = new TimeSeriesGenerator(100, 2_500, null);
        System.out.println("Done generating." + (System.currentTimeMillis()-tic));

        double[][][][] allTrajectories = tsg.getAllTrajectories();
        for (int classIdx = 0; classIdx < allTrajectories.length; classIdx++) {
            System.out.println(TimeSeriesGenerator.system_names[classIdx]);
            for (int trajectoryIdx = 0; trajectoryIdx < allTrajectories[classIdx].length; trajectoryIdx++) {
                double[][] ts = allTrajectories[classIdx][trajectoryIdx];

                // add noise
                double[][] tsUncorrected = TimeSeriesGenerator.addNoiseUncorrected(ts, 0.05, 200);

                DRQA drqa = new DRQA(ts, PhaseSpaceDistribution.maxPhaseSpaceDiameter(ts, ts) * 0.05);
                double[] rqa = drqa.allRQAMeasures(DRQA.STANDARD_RQA);
                System.out.println(Arrays.toString(rqa));
            }
        }
    }

    /**
     * Create one, two, and three-dimensional noise. Write to file.
     * Using the fixed recurrence rate criterion, the RQA measures seem to be independent from the dimensionality of the noise.
     * Using the diamter criterion, all three trajectories expose strongly different RQA measures, suspectedly simply because
     *   different dimensionalities imply strongly different recurrence chances.
     *   The differences in RR are expected to propagate to other RQA measures.
     */
    public void testNoiseLevels(){

        int tsLength = 20000;
        TimeSeriesGenerator tsg = new TimeSeriesGenerator(6, tsLength, "noise");
        double[][][][] noiseComponentes = tsg.getAllTrajectories();

        double[][][] noiseTrajectories = new double[][][]{
                new double[][]{noiseComponentes[7][0][0]},
                new double[][]{noiseComponentes[7][1][0], noiseComponentes[7][2][0]},
                new double[][]{noiseComponentes[7][3][0], noiseComponentes[7][4][0], noiseComponentes[7][5][0]},
        };

        for (int i = 0; i < 3; i++) {
            String filename = String.format("data/noise-dimensionality/noise-%sd.txt", i+1);
            try{
                Files.write(Paths.get(filename), TimeSeriesGenerator.trajectoryToCsvString(noiseTrajectories[i]).toString().getBytes());
            } catch (IOException e) { e.printStackTrace(); }
        }


    }



}
