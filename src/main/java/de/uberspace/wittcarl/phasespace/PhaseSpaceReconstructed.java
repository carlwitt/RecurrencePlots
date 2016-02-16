package de.uberspace.wittcarl.phasespace;

import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

import java.io.*;

/**
 * Applies Time Delay Embedding to univariate trajectories and estimates the parameters using the TISEAN software package.
 * The binaries of the programs mutual and false_nearest need to be provided in the folder external_binary.
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 17.08.15.
 */
public class PhaseSpaceReconstructed {

    public int embedding_dimension = 1, embedding_delay = 1, min_embedding_delay = -1;
    public double[][] trajectory;

    /** Stop increasing the embedding dimension once the fraction of false nearest neighbors drops below fnn_threshold. */
    public static double fnn_threshold = 0.05;
    int min_dimension = 1;

    public PhaseSpaceReconstructed(double[] trajectory, int minDimension) {
        this.min_dimension = minDimension;
        try {
            this.trajectory = embed(trajectory, false);
        } catch (IOException e) {
            this.trajectory = new double[][]{};
            e.printStackTrace();
        }
    }

    /**
     * Reconstruct the phase space trajectory using partial parameters. Given parameters are used, the others are estimated.
     * @param values The time series.
     * @param dim if <= 0, use the FNN algorithm with standard threshold. otherwise use the given value as dimension.
     * @param tau if <= 0, use the mutual information method to determine the embedding delay. otherwise use the given value as delay.
     */
    public PhaseSpaceReconstructed(double[] values, int dim, int tau) {

        try {
            if(dim > 0 && tau > 0)
                this.trajectory = PhaseSpaceReconstructed.embed(values, dim, tau);
            else{
                // estimate embedding parameters
                this.trajectory = embed(values, false);
                // no dimension given: replace with estimate
                if(dim <= 0) dim = embedding_dimension;
                // no delay given: replace with estimate
                if(tau <= 0) tau = embedding_delay;
                // repeat embedding using the final parameters
                if(!(dim == embedding_dimension && tau == embedding_delay))
                    this.trajectory = PhaseSpaceReconstructed.embed(values, dim, tau);
            }
        } catch (IOException e) {
            this.trajectory = new double[][]{};
            e.printStackTrace();
        }

    }

    /**
     * Reconstruct a phase space trajectory from the given univariate time series using time delay embedding.
     *
     * @param values univariate time series
     * @param dim    embedding dimension
     * @param tau    embedding delay
     * @return embedded[i][j] refers to the i-th dimension of the j-th phase space trajectory point
     */
    public static double[][] embed(double[] values, int dim, int tau) {
        double[][] embedded = new double[dim][values.length - (dim - 1) * tau];
        for (int i = 0; i < embedded[0].length; i++)
            for (int j = 0; j < dim; j++)
                embedded[j][i] = values[i + j * tau];

        return embedded;
    }

    /**
     * @param trajectory the univariate time series to embed
     * @param smallestLocalMinimum this option should always be true, false is maintained for reproducibility
     * @return the embedded trajectory, according to the estimated dimension and delay
     * @throws IOException
     */
    protected double[][] embed(double[] trajectory, boolean smallestLocalMinimum) throws IOException {

        Runtime rt = Runtime.getRuntime();

        // the string representation of the trajectory is fed to the executables via piping (faster and simpler than writing to disk and passing the path)
        String pipeContent = TimeSeriesGenerator.trajectoryToCsvString(new double[][]{trajectory}).toString();

        // determine delay as the first local minimum of the mutual information function
        final String mi_command = String.format("./external_binary/mutual -V0 -D%s", trajectory.length/2);
        Process ps_mi = rt.exec(mi_command);

        // feed the program the time series information via a pipe
        BufferedWriter out_mi_writer = new BufferedWriter(new OutputStreamWriter(ps_mi.getOutputStream()));
        try {
            out_mi_writer.write(pipeContent);
        } catch (IOException e){
            e.printStackTrace();
            System.err.println("Error in calling mutual information binary.");
        }
        out_mi_writer.close();

        // read and parse the program output
        BufferedReader reader = new BufferedReader(new InputStreamReader(ps_mi.getInputStream()));
        String row;
        double prevMI = Double.POSITIVE_INFINITY;
        double prevprevMI = Double.POSITIVE_INFINITY;
        double currentLocalMinimum = Double.POSITIVE_INFINITY;

        while ((row = reader.readLine()) != null) {
            if (row.startsWith("#shannon=")) continue;
            final String[] rowValues = row.split(" ");
            int delay = Integer.parseInt(rowValues[0]);
            double mutual_information = Double.parseDouble(rowValues[1]);
            if (mutual_information > prevMI && prevMI < prevprevMI  // is local minimum
                    && prevMI < currentLocalMinimum ) { // is smallest local minimum
                embedding_delay = delay - 1;
                if(min_embedding_delay < 0) min_embedding_delay = delay -1; // record first local minimum
                if(!smallestLocalMinimum) break;
                currentLocalMinimum = prevMI;
            }
            prevprevMI = prevMI;
            prevMI = mutual_information;
        }
        reader.close();

        // determine embedding dimension using the false nearest neighbors method
        int theiler_corrector = 1;
        int max_dimension = 10;
        final String fnn_command = String.format("./external_binary/false_nearest -M1,%s -d%s -t%s -V0", max_dimension, smallestLocalMinimum ? min_embedding_delay : embedding_delay, theiler_corrector);
        Process ps_fnn = rt.exec(fnn_command);
        BufferedWriter out_fnn_writer = new BufferedWriter(new OutputStreamWriter(ps_fnn.getOutputStream()));
        try {
            out_fnn_writer.write(pipeContent);
        } catch (IOException e){ System.err.println("Error in calling false_nearest binary."); }
        out_fnn_writer.close();

        embedding_dimension = min_dimension;
        reader = new BufferedReader(new InputStreamReader(ps_fnn.getInputStream()));
        while ((row = reader.readLine()) != null) {
            final String[] rowValues = row.split(" ");

            int dim = Integer.parseInt(rowValues[0]);
            double amount_fnn = Double.parseDouble(rowValues[1]);

            if (amount_fnn < fnn_threshold) {
                embedding_dimension = Math.max(min_dimension, dim);
                reader.close();
                break;
            }
        }

        return embed(trajectory, embedding_dimension, embedding_delay);

    }

}
