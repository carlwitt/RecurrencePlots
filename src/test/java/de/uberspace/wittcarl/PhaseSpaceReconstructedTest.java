package de.uberspace.wittcarl;

import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import de.uberspace.wittcarl.executable.BaseExecutable;
import de.uberspace.wittcarl.phasespace.PhaseSpaceReconstructed;
import junit.framework.TestCase;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 17.08.15.
 */
public class PhaseSpaceReconstructedTest extends TestCase{

    public void testReconstruct(){

        // read original trajectory
        final String file = "data/junit-test-data/lorenz-noisy-10000/lorenz-10000-0.0.csv";
        double[][] trajectory = BaseExecutable.readFile(file, ",");

        // reconstruct the phase space trajectory from first component and write the result to disk
        PhaseSpaceReconstructed reconstructed = new PhaseSpaceReconstructed(trajectory[0], 1);
        String reco = TimeSeriesGenerator.trajectoryToCsvString(reconstructed.trajectory).toString();
        BaseExecutable.writeStringToFile("data/junit-test-data/lorenz-10k-reconstructed-tde.csv", reco);

        // print RQA measures of original trajectory
        DRQA drqa_original = new DRQA(trajectory, trajectory, 0.05);
        drqa_original.computeRQA(2, 2, 2);

        // print RQA measures of reconstructed trajectory
        DRQA drqa = new DRQA(reconstructed.trajectory, reconstructed.trajectory, 0.05);
        drqa.computeRQA(2, 2, 2);

        drqa_original.printComparison(drqa);

    }


}