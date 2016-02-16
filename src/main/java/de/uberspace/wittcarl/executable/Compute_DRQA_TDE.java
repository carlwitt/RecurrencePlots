package de.uberspace.wittcarl.executable;

/**
 * Computes the definite recurrence quantification measures (DRQA) [1] for a given time series.
 * The time series is expected to be stored in text format, one observation per line.
 * If the time series is multivariate, use tabs as delimiters and pass the dimensionality as embedding dimension. No embedding will be applied then.
 *
 * An overview of the indices used:
 *  The points on the phase space trajectory are assigned zero based indices.
 *  The recurrence matrix at row i and column j refers to the distance between the i-th and j-th point on the phase space trajectory.
 *
 * [1] Carl Witt, "Recurrence Plot Clustering", Masters Thesis, Humboldt-Universiät zu Berlin, 2015.
 *
 */

import com.beust.jcommander.ParameterException;
import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.phasespace.PhaseSpaceReconstructed;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.JCommander;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Computes the Recurrence Plot and its RQA Measures for a one-dimensional time series by either applying the given time delay embedding parameters or estimating them.
 * input_dir: all files in this directory will be used as input
 * file_prefix: only those with a given prefix will be passed to the computation
 * threshold_fraction_pstd: specifies the recurrence threshold as a fraction of the maximum phase space trajectory diameter
 * embedding_dimension:
 *      if given an integer <= 0, will use the FNN algorithm with standard threshold to determine the embedding dimension.
 *      if given an integer  > 0, will use the given value as dimension.
 *      if given a double   >= 0, will use the FNN algorithm with the given threshold to determine the embedding dimension.
 * embedding_delay:
 *      if given an integer <= 0, will use the mutual information method to determine the embedding delay.
 *      if given an integer  > 0, will use the given value as delay.
 * write_rp: any value will be interpreted as true. outputs the recurrence plot as png image. this will require a lot of memory for long time series.
 * cross_recurrence_plot: any value will be interpreted as true. the input is now expected to have two columns (first two used).
 *      both columns are interpreted as 1-dimensional trajectories, are embedded and then the cross recurrence plot is computed.
 *      if the embedding parameters are to be estimated, they will be averaged to get phase space trajectories of the same length and dimensionality.
 * crt_limit: the maximal sum of two white lines to record in the CRT histogram.
 */
public class Compute_DRQA_TDE extends BaseExecutable {

    @Parameter(names = {"-i", "--in"}, description = "Name of the input directory.")
    String input_dir = "./";

    @Parameter(names = {"-f", "--prefix"}, description = "All files in the input directory with the given prefix will be matched.")
    String prefix = "";

    @Parameter(names = {"-p", "--out"}, description = "Output file for computed RQA measures.")
    String outputfile = "results.tsv";

    @Parameter(names = {"-d", "--out-dir"}, description = "Output directory for Recurrence Plots.")
    String out_dir = "./";

    @Parameter(names = {"-e", "--thresh-rel"}, description = "Recurrence threshold, relative to the maximum phase space diameter.")
    double fracMaxPSDDiameter = 0.05;

    @Parameter(names = {"-m", "--dim"}, description = "Embedding dimension. " +
            "An integer > 0 will be used as dimension. " +
            "If given an integer <= 0, the False Nearest Neighbors (FNN) algorithm with 5% threshold will be used to determine the dimension. " +
            "If given a fractional number >= 0, the FNN algorithm with the given threshold will be used to determine the embedding dimension.")
    double dimArgument = 1;

    @Parameter(names = {"-t", "--delay"}, description = "Embedding delay. If given an integer <= 0, the first minimum of the mutual information will be used.")
    int tau = 1;

    @Parameter(names = {"-w", "--write-rp"}, description = "Output Recurrence Plots as png files. Slow. Requires lots of memory for long time series.")
    boolean writeRP = false;

    @Parameter(names = {"-c", "--crp"}, description = "Compute the Cross Recurrence Plot instead of the Recurrence Plot. In this case, the two first two columns are interpreted as 1D input time series.")
    boolean crossRecurrencePlot = false;

    @Parameter(names = "--crt-limit", description = "The size of the Conditional Recurrence Time (CRT) histogram. If <= 1, no CRT is computed.")
    int crt_limit = 1;

    @Parameter(names = "--write-crt", description = "Plot the Conditional Recurrence Time histogram and save it as a png file.")
    boolean writeCRTPlot = false;

    @Parameter(names = {"-h", "--threads"}, description = "The number of threads to use. Parallelization is taking place on the file level, i.e. processing a single file is not sped up.")
    int num_threads = Runtime.getRuntime().availableProcessors()-1;

    public static void main(String[] args) throws IOException {

        final Compute_DRQA_TDE executable = new Compute_DRQA_TDE();

        try{
            JCommander argumentParser = new JCommander(executable, args);
            if(args.length == 0){
                argumentParser.setProgramName("DRQA");
                argumentParser.usage();
                return;
            }

        } catch (ParameterException e){
            e.printStackTrace();
        }


        DRQA.conditional_ww_limit = executable.crt_limit;
        int num_threads = executable.num_threads;

        // integer given, interpret as dimension
        final int dimension = executable.dimArgument % 1 == 0 ? (int) executable.dimArgument : 0;
        if(executable.dimArgument % 1 != 0) PhaseSpaceReconstructed.fnn_threshold = executable.dimArgument;   // double given, interpret as FNN threshold

        System.out.println((executable.crossRecurrencePlot ? "CROSS RECURRENCE PLOT" : "RECURRENCE PLOT"));
        System.out.println(String.format("Input Directory: %s", executable.input_dir));
        System.out.println(String.format("File Prefix: %s", executable.prefix));
        System.out.println(String.format("Recurrence Threshold: %s of max. phase space diameter", executable.fracMaxPSDDiameter*100));
        System.out.println(String.format("Dimension: %s", dimension == 0 ? String.format("Estimate using FNN < %.1f%%", PhaseSpaceReconstructed.fnn_threshold * 100) : dimension));
        System.out.println(String.format("Delay: %s", executable.tau == 0 ? "Estimate using mutual information." : executable.tau));


        System.out.println("Creating Tasks...");

        List<Callable<String>> tasks = new ArrayList<Callable<String>>();
        List<File> filteredFiles = getFilteredFiles(executable.input_dir, executable.prefix);
        for (final File file : filteredFiles) {
            Callable<String> c = new Callable<String>() { @Override public String call() {
                try { return processFile(file, executable.out_dir, dimension, executable.tau, executable.fracMaxPSDDiameter, executable.writeRP, executable.crossRecurrencePlot, executable.writeCRTPlot); }
                catch (IOException e) { e.printStackTrace(); return null; }
            } };
            tasks.add(c);
        }

        ExecutorService exec = Executors.newFixedThreadPool(num_threads);
        System.out.println(String.format("Processing %s tasks in %s threads.", tasks.size(), num_threads));

        BufferedWriter writer = Files.newBufferedWriter(Paths.get(executable.outputfile), Charset.defaultCharset(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
        writer.write(DRQA.separatedHeader()+"\n");

        try {
            long started = System.currentTimeMillis();
            List<Future<String>> results = exec.invokeAll(tasks);
            int written_lines = 0;
            for (Future<String> fr : results) {
                String line = fr.get();
                if (line == null) continue;
                writer.write(fr.get());
                writer.write("\n");
                written_lines++;
            }
            long elapsed_seconds = (System.currentTimeMillis() - started)/1000;
            System.out.println(String.format("Elapsed time: %sh:%sm", elapsed_seconds/60/60, (elapsed_seconds/60)%60));
            System.out.println(String.format("%s line(s) written to %s", written_lines, executable.outputfile));
        } catch (Exception e) {
            e.printStackTrace();
        }
        finally {
            exec.shutdown();
            writer.close();
        }


    }

    private static String processFile(File file, String out_dir, int dim, int tau, double fracMaxPSDDiameter, boolean writeRP, boolean crossRecurrencePlot, boolean writeCRTPlot) throws IOException {

        final String absolutePath = file.getAbsolutePath();
        final String filename = file.getName();

        System.out.println(String.format("Processing: %s", absolutePath));

        double[][] columns = BaseExecutable.readFile(absolutePath, ",");
        // the readfile method gives an empty array if the file is not a text file or something else goes wrong.
        if(columns.length == 0) return null;

        double[] values1, values2;
        double[][] embedded1, embedded2;
        PhaseSpaceReconstructed reconstructed;

        if(!crossRecurrencePlot){

            if(columns.length > 1) System.err.println(String.format("Using only first column of input file. %s discarded.", columns.length-1));
            values1 = columns[0];
            reconstructed = new PhaseSpaceReconstructed(values1, dim, tau);
            embedded1 = reconstructed.trajectory;
            embedded2 = reconstructed.trajectory;
            int dim1 = reconstructed.embedding_dimension;
            int tau1 = reconstructed.embedding_delay;

            if(dim <= 0) System.out.println(String.format("Estimated embedding dimension: %s", dim1));
            if(tau <= 0) System.out.println(String.format("Estimated time delay: %s", tau1));

            // finally used parameters
            dim = dim <= 0 ? dim1 : dim;
            tau = tau <= 0 ? tau1 : tau;

        } else {
            if(columns.length > 2) System.err.println(String.format("Using only first two columns of input file. %s discarded.", columns.length-2));
            values1 = columns[0];
            values2 = columns[1];

            // use time delay embedding and possibly estimation of dim and tau (if values <= 0 are given)
            reconstructed = new PhaseSpaceReconstructed(values1, dim, tau);
            int dim1 = dim <= 0 ? reconstructed.embedding_dimension : dim;
            int tau1 = tau <= 0 ? reconstructed.embedding_delay : tau;
            PhaseSpaceReconstructed reconstructed2 = new PhaseSpaceReconstructed(values2, dim, tau);
            int dim2 = dim <= 0 ? reconstructed2.embedding_dimension : dim;
            int tau2 = tau <= 0 ? reconstructed2.embedding_delay : tau;

            // if embedding parameters differ, repeat embedding with the compromise parameters
            int avg_dim=dim1, avg_tau=tau1;
            if(dim1 != dim2) avg_dim = Math.max(1, (dim1+dim2)/2);
            if(tau1 != tau2) avg_tau = Math.max(1, (tau1+tau2)/2);
            if(avg_dim != dim1 || avg_tau != tau1) reconstructed = new PhaseSpaceReconstructed(values1, avg_dim, avg_tau);
            if(avg_dim != dim2 || avg_tau != tau2) reconstructed2 = new PhaseSpaceReconstructed(values2, avg_dim, avg_tau);

            // if estimation was required, print the estimation results
            if(dim <= 0) System.out.println(String.format("Estimated embedding dimension: %s (first column) %s (second column) %s (used for both)", dim1, dim2, avg_dim));
            if(tau <= 0) System.out.println(String.format("Estimated embedding tau: %s (first column) %s (second column) %s (used for both)", tau1, tau2, avg_tau));

            embedded1 = reconstructed.trajectory;
            embedded2 = reconstructed2.trajectory;

            // finally used parameters
            dim = dim <= 0 ? avg_dim : dim;
            tau = tau <= 0 ? avg_tau : tau;
        }



        // record start time
        long started = System.currentTimeMillis();

        DRQA drqa = new DRQA(embedded1, embedded2, fracMaxPSDDiameter);
        System.out.println(String.format("DRQA computed in: %s ms", System.currentTimeMillis() - started));

        drqa.computeRQA(2, 2, 2);

        String line = drqa.toSeparatedString(filename, TimeSeriesGenerator.findClassName(filename), DRQA.EmbeddingMethod.ORIGINAL_TRAJECTORY.ordinal(), dim, tau, DRQA.LMinMethod.L_MIN_FIX.ordinal(), 2, 0, embedded1[0].length);

//        Files.write(Paths.get(filename+"_CRT.tsv"), drqa.conditionalRecurrenceTsv().toString().getBytes());
        if(writeCRTPlot) drqa.writeCRTImages(String.format("%s/%s_%s_dim%stau%s_signature.png", out_dir, filename, crossRecurrencePlot?"CRP":"RP", dim, tau));

        if(writeRP) drqa.writeRP(embedded1, embedded2, fracMaxPSDDiameter, String.format("%s/%s_%s_dim%stau%s.png", out_dir, filename, crossRecurrencePlot?"CRP":"RP", dim, tau));

        return line;
    }

    /**
     * The idea was to lower bound the number of steps before the trajectory can return to a given point.
     * When computing all recurrences for single point (computing the recurrence matrix column-wise), the in-between cells can be skipped.
     * But the overhead is slowing things down, since the computation of the single cells are not really expensive. Unfolding of loops is also not possible anymore.
     * For the Lorenz System, the following values emerged:
     *      maxAdjDistance: 12.190735488951919
     *      average adjacent distance: 4.725239546436434
     *      maxSpaceDistance: 41.37204034066555
     * So the saving was even for the farthest points on the trajectory only four cells to skip.
     * I tried the following:
     *      earliestPossibleRecurrenceIn = Math.max(1, (int) Math.ceil(distance / maxAdjDistance));
     *      column_index += earliestPossibleRecurrenceIn;
     */
    private static void maxAdjacentDistance(int dim, double[][] embedded) {
        // compute maximum phase space distance between any two subsequent points on the trajectory
        double ssoq;
        double maxAdjDistance = 0;
        double sumAdjDistance = 0;
        for (int i = 1; i < embedded.length; i++) {
            // compute euclidean distance in variable number of dimensions
            ssoq = 0;
            for (int k = 0; k < dim; k++) {
                ssoq += (embedded[i][k] - embedded[i-1][k]) * (embedded[i][k] - embedded[i-1][k]);
            }
            final double distance = Math.sqrt(ssoq);
            sumAdjDistance += distance;
            if(distance > maxAdjDistance) maxAdjDistance = distance;
        }
        System.out.println(String.format("maxAdjDistance: %s", maxAdjDistance));
        final double avgAdjDistance = sumAdjDistance / (embedded.length - 1);
        System.out.println(String.format("average adjacent distance: %s", avgAdjDistance));
    }
}
