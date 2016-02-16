package de.uberspace.wittcarl;

import de.uberspace.wittcarl.datagenerator.FlowSystem;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import junit.framework.TestCase;
import org.junit.Assert;

import java.util.Arrays;

public class FlowSystemTest extends TestCase {

    // check whether the interpolation of the continuous output model works as expected.
    public void testInterpolation() throws Exception {

        FlowSystem flowSystem = new FlowSystem(TimeSeriesGenerator.LorenzEquations.standard(), new double[]{1,2,3}, 0.005, 5);

        // gives the original phase space trajectory
        final double[][] base_trajectory = flowSystem.getTrajectory(0.005);
        System.out.println(String.format("flowSystem.getTransformedTrajectory(0.005):\n%s", TimeSeriesGenerator.trajectoryToCsvString(base_trajectory)));

        // gives the original phase space trajectory plus linear interpolations in between.
        final double[][] interpolated = flowSystem.getTrajectory(0.0025);
        System.out.println(String.format("flowSystem.getTransformedTrajectory(0.0025):\neach second step should be equal to above output. The values in between should be linearly interpolated\n%s", TimeSeriesGenerator.trajectoryToCsvString(interpolated)));

        for (int i = 0; i < interpolated.length; i++) {
            assertEquals(base_trajectory[0][i], interpolated[0][i*2]);
            assertEquals(base_trajectory[1][i], interpolated[1][i*2]);
            assertEquals(base_trajectory[2][i], interpolated[2][i*2]);
        }
    }

    /**
     * Assert that cutTrajectories yields the same as manually cutting the trajectories.
     */
    public void testCutTrajectories(){

        FlowSystem roessler = new FlowSystem(TimeSeriesGenerator.RoesslerEquations.standard(), new double[]{1,2,3}, 0.001, 11);
        final double[][] base_trajectory = roessler.getTrajectory();

        double[][] expectedFirst = new double[][]{Arrays.copyOfRange(base_trajectory[0],1,6), Arrays.copyOfRange(base_trajectory[1],1,6), Arrays.copyOfRange(base_trajectory[2],1,6)};
        double[][] expectedSecnd = new double[][]{Arrays.copyOfRange(base_trajectory[0],6,11), Arrays.copyOfRange(base_trajectory[1],6,11), Arrays.copyOfRange(base_trajectory[2],6,11)};

        double[][][] trajectories = FlowSystem.cutTrajectories(TimeSeriesGenerator.RoesslerEquations.standard(), new double[]{1,2,3}, 0.001, 2, 5, 1);

        Assert.assertArrayEquals(expectedFirst, trajectories[0]);
        Assert.assertArrayEquals(expectedSecnd, trajectories[1]);

    }

}