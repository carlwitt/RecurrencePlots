package de.uberspace.wittcarl.datagenerator;

import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 13.08.15.
 */


/**
 * The base class for the Rössler and Lorenz Systems.
 * Defines the integration algorithm, fourth order Runge-Kutta (standard in this domain).
 * Stores the trajectory data and handles interpolation of the trajectory data.
 */
public class FlowSystem {

    int dimensionality;
    ContinuousOutputModel continuousOutputModel = new ContinuousOutputModel();
    double originalStepsize;

    public FlowSystem(FirstOrderDifferentialEquations equations, double[] initialConditions, double stepSize, int numSteps) {

        dimensionality = equations.getDimension();
        originalStepsize = stepSize;

        // using fourth-order Runge-Kutta numerical integration
        ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(stepSize);
        integrator.addStepHandler(continuousOutputModel);

        // run integrator
        integrator.integrate(equations,
                0,                              // start time
                initialConditions,              // initial conditions
                numSteps*stepSize,              // end time
                new double[dimensionality]);    // storage

    }

    /**
     * @return all computed data, without any interpolation.
     */
    public double[][] getTrajectory(){ return getTrajectory(originalStepsize); }

    /**
     * @param stepSize ∆t
     * @return The data, interpolated to evenly spaced observations at i*∆t
     */
    public double[][] getTrajectory(double stepSize){
        int length = (int) Math.floor(continuousOutputModel.getFinalTime() / stepSize);
        double[][] interpolated_trajectory = new double[dimensionality][length];
        for (int i = 0; i < length; i++) {
            continuousOutputModel.setInterpolatedTime(i * stepSize);
            final double[] interpolatedState = continuousOutputModel.getInterpolatedState();
            for (int dim = 0; dim < dimensionality; dim++)
                interpolated_trajectory[dim][i] = interpolatedState[dim];
        }
        return interpolated_trajectory;
    }

    public double[][] getTrajectory(int settling_time, double stepSize, int chunk_idx, int timeSeriesLength){
        double fromTime = chunk_idx*timeSeriesLength*stepSize;
        double toTime = (chunk_idx+1)*timeSeriesLength*stepSize;
        double[][] interpolated_trajectory = new double[dimensionality][timeSeriesLength];
        for (int i = 0; i < timeSeriesLength; i++) {
            continuousOutputModel.setInterpolatedTime(settling_time*stepSize + fromTime + i*stepSize);
            final double[] interpolatedState = continuousOutputModel.getInterpolatedState();
            for (int dim = 0; dim < dimensionality; dim++)
                interpolated_trajectory[dim][i] = interpolatedState[dim];
        }
        return interpolated_trajectory;
    }

    public static double[][][] cutTrajectories(FirstOrderDifferentialEquations equations, double[] initialConditions, double stepSize, int numberOfTimeSeries, int timeSeriesLength, int settling_time) {
        FlowSystem flowSystem = new FlowSystem(equations, initialConditions, stepSize, numberOfTimeSeries * timeSeriesLength + settling_time);
        double[][][] trajectories = new double[numberOfTimeSeries][equations.getDimension()][timeSeriesLength];
        for (int i = 0; i < numberOfTimeSeries; i++) {
            trajectories[i] = flowSystem.getTrajectory(settling_time, stepSize, i, timeSeriesLength);
        }
        return trajectories;
    }

}
