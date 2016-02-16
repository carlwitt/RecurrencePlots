package de.uberspace.wittcarl.experimental;

import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import de.uberspace.wittcarl.executable.BaseExecutable;
import de.uberspace.wittcarl.phasespace.PhaseSpaceReconstructed;
import junit.framework.TestCase;

import java.io.IOException;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 20.08.15.
 */
public class NoiseAndEmbeddingTest extends TestCase {

    /**
     * The problem that for noise an embedding dimension of one is estimated is also present for 1% and 0.5% FNN threshold (although the dimensionalities without noise increase)
     * @throws IOException
     */
    public void testFNNThreshold() throws IOException {

        // generate a base trajectory of the standard lorenz system
        final int numSteps = 25000;
        PhaseSpaceReconstructed.fnn_threshold = 0.05;

        TimeSeriesGenerator timeSeriesGenerator = new TimeSeriesGenerator(3, numSteps, null);

        double[][][][] allTrajectories = timeSeriesGenerator.getAllTrajectories();

        for (int i = 0; i < allTrajectories.length-4; i++) {
            String name = TimeSeriesGenerator.system_names[i];
            System.out.println("\n"+name+"\n");

            for (double noiseRatio : new double[]{0.0, 0.01, 0.1}) {
                System.out.println(noiseRatio);

                for (int j = 0; j < allTrajectories[i].length; j++) {
                    double[][] original = allTrajectories[i][j];

                    double[][] noisyTrajectory = TimeSeriesGenerator.addNoise(original, 0.05, noiseRatio);
//                    DRQA drqa = new DRQA(noisyTrajectory, noisyTrajectory, 0.05);
                    PhaseSpaceReconstructed phaseSpaceReconstructed = new PhaseSpaceReconstructed(noisyTrajectory[0], 1);
                    System.out.println(String.format("phaseSpaceReconstructed.embedding_dimension: %s", phaseSpaceReconstructed.embedding_dimension));
                    System.out.println(String.format("phaseSpaceReconstructed.embedding_delay: %s", phaseSpaceReconstructed.embedding_delay));
                }

            }
        }
    }

    public void testDelayEstimationRoessler(){

//        double[][] trajectory = BaseExecutable.readFile("data/junit-test-data/pst_1_roessler_standard28.txt", ",");
        double[][] trajectory = BaseExecutable.readFile("data/junit-test-data/lorenz-noisy-10000/lorenz-10000-0.8.csv", ",");

        PhaseSpaceReconstructed phaseSpaceReconstructed = new PhaseSpaceReconstructed(trajectory[0], 1);
        System.out.println(phaseSpaceReconstructed.embedding_dimension);
        System.out.println(phaseSpaceReconstructed.embedding_delay);

    }

    // see how RQA measures change when switching from the original lorenz to its one dimensional projection
    public void testProject(){

        // read original trajectory
        final String file = "data/junit-test-data/lorenz-noisy-10000/lorenz-10000-0.0.csv";
        double[][] trajectory = BaseExecutable.readFile(file, ",");

        // project to first dimension
        double[][] projection = new double[][]{trajectory[0]};
        BaseExecutable.writeStringToFile("data/junit-test-data/lorenz-10k-x-component.csv", TimeSeriesGenerator.trajectoryToCsvString(projection).toString());

        // print RQA measures of original trajectory
        DRQA drqa_original = new DRQA(trajectory, trajectory, 0.05);
        drqa_original.computeRQA(2, 2, 2);

        // print RQA measures of projected trajectory
        DRQA drqa = new DRQA(projection, projection, 0.05);
        drqa.computeRQA(2, 2, 2);

        drqa_original.printComparison(drqa);

    }

}
