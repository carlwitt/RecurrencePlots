package de.uberspace.wittcarl.datagenerator;

import de.uberspace.wittcarl.phasespace.PhaseSpaceDistribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.Well19937a;
import org.apache.commons.math3.special.Gamma;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;
import java.util.logging.Logger;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Data: 03.08.15.
 * <p>
 * Time series generator for the evaluation framework.
 * For a given number of time series and observations, creates a data set consisting of nine classes, as described in Section 5.1.1 of the thesis.
 */
public class TimeSeriesGenerator {

    public static String[] system_names = new String[]{"oscillator", "roessler_standard", "roessler_funnel", "lorenz_standard", "lorenz_noisy_oscillation", "autoregressive_negative", "autoregressive_positive", "noise_normal", "noise_gamma"};
    public double[][][] lorenz_standard, lorenz_noisy_oscillation, roessler_standard, roessler_funnel;
    // first dimension denotes trajectory index, second dimension refers to dimension, third to time step
    public double[][][] oscillator;
    public double[][][] noise_normal, noise_gamma, autoregressive_negative, autoregressive_positive;

    public static void main(String[] args) throws IOException {

        if (args.length < 2) {
            System.out.println("Usage: numberOfTs lengthOfTs outDir [system] [noise_ratio]");
            System.exit(0);
        }

        int numberOfTS = Integer.parseInt(args[0]);
        int lengthOfTS = Integer.parseInt(args[1]);
        String outDir = args.length > 2 ? args[2] : "./";
        String system = args.length > 3 ? args[3] : null;
        double noiseRatio = args.length > 4 ? Double.parseDouble(args[4]) : 0.0;

        TimeSeriesGenerator tsg = new TimeSeriesGenerator(numberOfTS, lengthOfTS, system);
        tsg.writeOut(outDir, system);

    }

    public TimeSeriesGenerator(int timeSeriesLength) {
        this(1, timeSeriesLength, null);
    }

    /**
     * @param system one of oscillator, roessler, lorenz, autoregressive, noise or null for all
     */
    public TimeSeriesGenerator(int numberOfTimeSeries, int timeSeriesLength, String system) {

        System.out.println("Computing Data");
        if (system == null || system.equals("oscillator")) {
            System.out.println("Oscillator...");
            oscillator = Oscillator.sample(numberOfTimeSeries, timeSeriesLength);
        }

        final int settling_steps = 1000;
        if (system == null || system.equals("roessler")) {
            System.out.println("Rössler System...");
            final double dt = RoesslerEquations.DEFAULT_STEP_SIZE;
            final double[] initialConditions = RoesslerEquations.DEFAULT_INITIAL_CONDITIONS;
            roessler_standard = FlowSystem.cutTrajectories(RoesslerEquations.standard(), initialConditions, dt, numberOfTimeSeries, timeSeriesLength, settling_steps);
            roessler_funnel = FlowSystem.cutTrajectories(RoesslerEquations.funnel(), initialConditions, dt, numberOfTimeSeries, timeSeriesLength, settling_steps);
        }
        if (system == null || system.equals("lorenz")) {
            System.out.println("Lorenz System...");
            lorenz_standard = FlowSystem.cutTrajectories(LorenzEquations.standard(), LorenzEquations.DEFAULT_INITIAL_CONDITIONS, LorenzEquations.DEFAULT_STEP_SIZE, numberOfTimeSeries, timeSeriesLength, settling_steps);
            lorenz_noisy_oscillation = FlowSystem.cutTrajectories(LorenzEquations.noisyOscillation(), LorenzEquations.DEFAULT_INITIAL_CONDITIONS, LorenzEquations.DEFAULT_STEP_SIZE, numberOfTimeSeries, timeSeriesLength, settling_steps);
        }

        if (system == null || system.equals("autoregressive")) {
            System.out.println("Autoregressive...");
            autoregressive_negative = Autoregressive.sample(-0.95, numberOfTimeSeries, timeSeriesLength);
            autoregressive_positive = Autoregressive.sample(0.95, numberOfTimeSeries, timeSeriesLength);
        }
        if (system == null || system.equals("noise")) {
            System.out.println("Noise...");
            noise_normal = Noise.normal(numberOfTimeSeries, timeSeriesLength);
            noise_gamma = Noise.gamma(numberOfTimeSeries, timeSeriesLength);
        }

    }

    /**
     * Creates a copy of the trajectory with a normally distributed displacement added to each component.
     *
     * @param trajectory                     The phase space trajectory, first dimension is referring to the phase space dimension, second dimension to time step.
<<<<<<< HEAD:src/main/java/de/uberspace/wittcarl/datagenerator/TimeSeriesGenerator.java
     * @param referenceThresholdFractionPSTD Specifies the recurrence threshold, as fraction of of the maximum phase space diameter.
     * @param noiseRatio                     The ratio between expected distance from original point and recurrence threshold. A ratio of 0 corresponds to no noise. A ratio of 1 corresponds to an expected displacement equal to the recurrence threshold.
     * @return A copy of the input data, with noise added.
     */
    public static double[][] addNoise(double[][] trajectory, double referenceThresholdFractionPSTD, double noiseRatio) {
        final int dim = trajectory.length;

        // the reference threshold is a fraction (usually 5%) of the phase space diameter (maximum distance between any two points on the phase space trajectory)
        double trajectoryMaxDistance = PhaseSpaceDistribution.maxPhaseSpaceDiameter(trajectory, trajectory);
        double referenceThreshold = referenceThresholdFractionPSTD * trajectoryMaxDistance;

        // to assure that the expected distance between displaced point and original point equals a given value, use the expected value of the scaled chi distribution with dim degrees of freedom
        double standardDeviation = Math.sqrt(2) * noiseRatio * referenceThreshold * Gamma.gamma(0.5 * dim) / Gamma.gamma(0.5 * (dim + 1));
        return addNoiseWithSD(trajectory, noiseRatio, standardDeviation);
    }

    /**
     * Creates a copy of the trajectory with a normally distributed displacement added to each component. Doesn't normalize the magnitude of the dimension vectors according to dimension (as addNoise does).
     * @param trajectory                     The phase space trajectory, first dimension is referring to the phase space dimension, second dimension to time step.
     * @param referenceThresholdFractionPSTD Specifies the recurrence threshold, as fraction of of the maximum phase space diameter.
     * @param noiseRatio                     The ratio between expected distance from original point and recurrence threshold. A ratio of 0 corresponds to no noise. A ratio of 1 corresponds to an expected displacement equal to the recurrence threshold.
     * @return A copy of the input data, with noise added.
     */
    public static double[][] addNoiseUncorrected(double[][] trajectory, double referenceThresholdFractionPSTD, double noiseRatio) {
        // the reference threshold is a fraction (usually 5%) of the phase space diameter (maximum distance between any two points on the phase space trajectory)
        double trajectoryMaxDistance = PhaseSpaceDistribution.maxPhaseSpaceDiameter(trajectory, trajectory);
        double referenceThreshold = referenceThresholdFractionPSTD * trajectoryMaxDistance;

        double standardDeviation = noiseRatio * referenceThreshold;
        return addNoiseWithSD(trajectory, noiseRatio, standardDeviation);
    }

    /**
     * Creates a copy of the trajectory with a normally distributed displacement added to each component.
     * @param trajectory                     The phase space trajectory, first dimension is referring to the phase space dimension, second dimension to time step.
     * @param noiseRatio                     The ratio between expected distance from original point and recurrence threshold. A ratio of 0 corresponds to no noise. A ratio of 1 corresponds to an expected displacement equal to the recurrence threshold.
     * @param standardDeviation              The standard deviation of the normal distribution from which the components of the perturbation vector are drawn.
     * @return A copy of the input data, with noise added.
     */
=======
     * @param referenceThresholdFractionPSTD Specifies the recurrence threshold, as fraction of of the maximum phase space diameter.
     * @param noiseRatio                     The ratio between expected distance from original point and recurrence threshold. A ratio of 0 corresponds to no noise. A ratio of 1 corresponds to an expected displacement equal to the recurrence threshold.
     * @return A copy of the input data, with noise added.
     */
    public static double[][] addNoise(double[][] trajectory, double referenceThresholdFractionPSTD, double noiseRatio) {
        final int dim = trajectory.length;

        // the reference threshold is a fraction (usually 5%) of the phase space diameter (maximum distance between any two points on the phase space trajectory)
        double trajectoryMaxDistance = DRQA.maxPhaseSpaceDiameter(trajectory, trajectory);
        double referenceThreshold = referenceThresholdFractionPSTD * trajectoryMaxDistance;

        // to assure that the expected distance between displaced point and original point equals a given value, use the expected value of the scaled chi distribution with dim degrees of freedom
        double standardDeviation = Math.sqrt(2) * noiseRatio * referenceThreshold * Gamma.gamma(0.5 * dim) / Gamma.gamma(0.5 * (dim + 1));
        return addNoiseWithSD(trajectory, noiseRatio, standardDeviation);
    }

    /**
     * Creates a copy of the trajectory with a normally distributed displacement added to each component. Doesn't normalize the magnitude of the dimension vectors according to dimension (as addNoise does).
     * @param trajectory                     The phase space trajectory, first dimension is referring to the phase space dimension, second dimension to time step.
     * @param referenceThresholdFractionPSTD Specifies the recurrence threshold, as fraction of of the maximum phase space diameter.
     * @param noiseRatio                     The ratio between expected distance from original point and recurrence threshold. A ratio of 0 corresponds to no noise. A ratio of 1 corresponds to an expected displacement equal to the recurrence threshold.
     * @return A copy of the input data, with noise added.
     */
    public static double[][] addNoiseUncorrected(double[][] trajectory, double referenceThresholdFractionPSTD, double noiseRatio) {
        // the reference threshold is a fraction (usually 5%) of the phase space diameter (maximum distance between any two points on the phase space trajectory)
        double trajectoryMaxDistance = DRQA.maxPhaseSpaceDiameter(trajectory, trajectory);
        double referenceThreshold = referenceThresholdFractionPSTD * trajectoryMaxDistance;

        double standardDeviation = noiseRatio * referenceThreshold;
        return addNoiseWithSD(trajectory, noiseRatio, standardDeviation);
    }

    /**
     * Creates a copy of the trajectory with a normally distributed displacement added to each component.
     * @param trajectory                     The phase space trajectory, first dimension is referring to the phase space dimension, second dimension to time step.
     * @param noiseRatio                     The ratio between expected distance from original point and recurrence threshold. A ratio of 0 corresponds to no noise. A ratio of 1 corresponds to an expected displacement equal to the recurrence threshold.
     * @param standardDeviation              The standard deviation of the normal distribution from which the components of the perturbation vector are drawn.
     * @return A copy of the input data, with noise added.
     */
>>>>>>> origin/indefinite-lines-via-position:src/de/uberspace/wittcarl/datagenerator/TimeSeriesGenerator.java
    public static double[][] addNoiseWithSD(double[][] trajectory, double noiseRatio, double standardDeviation) {

        final int dim = trajectory.length;
        final int steps = trajectory[0].length;
        double[][] noisy_copy = new double[dim][];

        // if no noise should be added, just return a copy of the input data
        if (noiseRatio == 0.0) {
            for (int dim_idx = 0; dim_idx < dim; dim_idx++)
                noisy_copy[dim_idx] = Arrays.copyOf(trajectory[dim_idx], steps);
            return noisy_copy;
        }

        NormalDistribution normalDistribution = new NormalDistribution(new Well19937a(0), 0, standardDeviation);

        // the displacement is drawn from a spherical normal distribution, so the covariances are 0.
        // accordingly, each component of the vector can be displaced independently.
        for (int dim_idx = 0; dim_idx < dim; dim_idx++) {
            // copy the original data to the output data
            noisy_copy[dim_idx] = Arrays.copyOf(trajectory[dim_idx], steps);
            // displace each component's value by a random "error" as given
            for (int step = 0; step < steps; step++) noisy_copy[dim_idx][step] += normalDistribution.sample();
        }
        return noisy_copy;
    }

    public static String findClassName(String name) {
        for (String system_name : system_names) {
            if (name.contains(system_name)) return system_name;
        }
        return "NOT_FOUND";
    }

    /**
     * @return [system][trajectory index][dimension][time step]
     */
    public double[][][][] getAllTrajectories() {
        return new double[][][][]{oscillator, roessler_standard, roessler_funnel, lorenz_standard, lorenz_noisy_oscillation, autoregressive_negative, autoregressive_positive, noise_normal, noise_gamma};
    }

    public static class Oscillator {

        public static int minSamplesPerPeriod = 32;     // defines the minimum wavelength of an oscillation
        public static int minPeriods = 3;               // defines the maximum wavelength of an oscillation, e.g. 3 periods = 2 repitions of the first period

        public static double[][][] sample(int numberOfTimeSeries, int timeSeriesLength) {

            double[][][] trajectories = new double[numberOfTimeSeries][2][timeSeriesLength];

            double maxSamplesPerPeriod = timeSeriesLength / minPeriods;

            for (int ts_idx = 0; ts_idx < numberOfTimeSeries; ts_idx++) {
                double samplesPerPeriod = 1. * (ts_idx + 1) / numberOfTimeSeries * (maxSamplesPerPeriod - minSamplesPerPeriod) + minSamplesPerPeriod;
                for (int step = 0; step < timeSeriesLength; step++) {
                    trajectories[ts_idx][0][step] = Math.sin(2.0 * Math.PI * step / samplesPerPeriod);
                    trajectories[ts_idx][1][step] = Math.cos(2.0 * Math.PI * step / samplesPerPeriod);
                }
            }
            return trajectories;
        }

    }

    public static class RoesslerEquations implements FirstOrderDifferentialEquations {
        public static double DEFAULT_STEP_SIZE = 0.085;
        public static double[] DEFAULT_INITIAL_CONDITIONS = new double[]{1, 2, 3};
        public static double DEFAULT_A = 0.2, DEFAULT_B = 0.2, DEFAULT_C = 5.7;
        public static double FUNNEL_A = 0.2925, FUNNEL_B = 0.1, FUNNEL_C = 8.5;
        public double a, b, c;

        public static RoesslerEquations standard() {
            return new RoesslerEquations(DEFAULT_A, DEFAULT_B, DEFAULT_C);
        }

        public static RoesslerEquations funnel() {
            return new RoesslerEquations(FUNNEL_A, FUNNEL_B, FUNNEL_C);
        }

        public RoesslerEquations(double a, double b, double c) {
            this.a = a;
            this.b = b;
            this.c = c;
        }

        public int getDimension() {
            return 3;
        }

        public void computeDerivatives(double t, double[] y, double[] yDot) {
            yDot[0] = -y[1] - y[2];           // dx/dt = -y-z
            yDot[1] = y[0] + a * y[1];          // dy/dt = x+ay
            yDot[2] = b + y[2] * (y[0] - c);      // dz/dt = b+z(x-c)
        }
    }

    public static class LorenzEquations implements FirstOrderDifferentialEquations {
        public static double DEFAULT_STEP_SIZE = 0.025;
        public static double[] DEFAULT_INITIAL_CONDITIONS = new double[]{-8, 8, 27};
        public static double DEFAULT_SIGMA = 10, DEFAULT_BETA = 8. / 3, DEFAULT_RHO = 28;
        public static double NOISY_OSCILLATION_SIGMA = 10, NOISY_OSCILLATION_BETA = 8. / 3, NOISY_OSCILLATION_RHO = 198;

        public double sigma, beta, rho;

        public static LorenzEquations standard() {
            return new LorenzEquations(DEFAULT_SIGMA, DEFAULT_BETA, DEFAULT_RHO);
        }

        public static LorenzEquations noisyOscillation() {
            return new LorenzEquations(NOISY_OSCILLATION_SIGMA, NOISY_OSCILLATION_BETA, NOISY_OSCILLATION_RHO);
        }

        public LorenzEquations(double sigma, double beta, double rho) {
            this.sigma = sigma;
            this.beta = beta;
            this.rho = rho;
        }

        public int getDimension() {
            return 3;
        }

        public void computeDerivatives(double t, double[] y, double[] yDot) {
            yDot[0] = sigma * (y[1] - y[0]);              // dx/dt = sigma(y-x)
            yDot[1] = y[0] * (rho - y[2]) - y[1];       // dy/dt = x(rho-z)-y
            yDot[2] = y[0] * y[1] - beta * y[2];        // dz/dt = xy - beta z
        }
    }

    public static class Autoregressive {
        public static double[][][] sample(double correlation, int numberOfTimeSeries, int timeSeriesLength) {
            double[][][] trajectories = new double[numberOfTimeSeries][1][timeSeriesLength];
            NormalDistribution normalDistribution = new NormalDistribution(0, 0.01);
            for (int ts_idx = 0; ts_idx < numberOfTimeSeries; ts_idx++) {
                double lastValue = normalDistribution.sample();
                for (int step = 0; step < timeSeriesLength; step++) {
                    trajectories[ts_idx][0][step] = correlation * lastValue + normalDistribution.sample();
                    lastValue = trajectories[ts_idx][0][step];
                } // for step
            } // for time series
            return trajectories;
        } // sample()
    } // Autoregressive

    public static class Noise {

        public static double[][][] normal(int numberOfTimeSeries, int timeSeriesLength) {
            return sample(new NormalDistribution(0, 1), numberOfTimeSeries, timeSeriesLength);
        }

        public static double[][][] gamma(int numberOfTimeSeries, int timeSeriesLength) {
            return sample(new GammaDistribution(0.5, 1), numberOfTimeSeries, timeSeriesLength);
        }

        static double[][][] sample(AbstractRealDistribution distribution, int numberOfTimeSeries, int timeSeriesLength) {
            double[][][] trajectories = new double[numberOfTimeSeries][1][timeSeriesLength];
            for (int i = 0; i < numberOfTimeSeries; i++)
                trajectories[i][0] = distribution.sample(timeSeriesLength);
            return trajectories;
        }
    }

    /**
     * Writes the text versions of each trajectory to a file each.
     *
     * @param outDir The prefix for each filename, must end with a valid directory separator, e.g. "out/"
     * @throws IOException If something goes wrong writing the text file.
     */
    public void writeOut(String outDir, String system) throws IOException {

        double[][][][] trajectories_by_system = getAllTrajectories();

        Logger.getGlobal().info("Writing time series to text files.");

        // how many digits for the index numbers in the file name (keeping them at the same number of digits improves lexicographic sorting)
        int num_digits_class = (int) Math.ceil(Math.log10(trajectories_by_system.length));
        int num_digits_timeseries = (int) Math.ceil(Math.log10(trajectories_by_system[0].length));

        // for each trajectories set
        for (int i = 0; i < trajectories_by_system.length; i++) {
            double[][][] trajectories = trajectories_by_system[i];
            String name = system_names[i];

            if (system != null && !name.startsWith(system)) continue;

            System.out.println(name);
            // for each trajectory
            for (int trajectory_idx = 0; trajectory_idx < trajectories.length; trajectory_idx++) {
                String filename = String.format(String.format("pst_%%0%sd_%%s%%0%dd.txt", num_digits_class, num_digits_timeseries), i, name, trajectory_idx);
                Files.write(Paths.get(outDir + filename), trajectoryToCsvString(trajectories[trajectory_idx]).toString().getBytes());
            }
        }

    }

    /**
     * Creates a string that contains one phase space trajectory point per line, the coordinates in the different dimensions separated by commas.
     *
     * @param trajectory The phase space trajectory.
     * @return The string builder containing the string.
     */
    public static StringBuilder trajectoryToCsvString(double[][] trajectory) {

        final int dimensions = trajectory.length;
        final int steps = trajectory[0].length;

        StringBuilder builder = new StringBuilder(trajectory.length * steps * 20);
        Locale.setDefault(Locale.US);
        for (int step = 0; step < steps; step++) {
            for (int dim = 0; dim < dimensions; dim++) {
                double[] dimension = trajectory[dim];
                builder.append(String.format("%.16f%s", dimension[step], dim < dimensions - 1 ? "," : ""));
            }
            builder.append("\n");
        }
        return builder;
    }
}
