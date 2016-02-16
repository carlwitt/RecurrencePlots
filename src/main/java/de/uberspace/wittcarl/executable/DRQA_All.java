package de.uberspace.wittcarl.executable;

import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.phasespace.PhaseSpaceReconstructed;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

import java.io.BufferedWriter;
import java.io.Console;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Takes a set of input files and computes the DRQA measures in parallel.
 * Simulates the effect of noise and embedding (treats each input file with each combination of noise level {@link DRQA_All#noiseRatiosById} and embedding method {@link de.uberspace.wittcarl.DRQA.EmbeddingMethod}
 *
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 17.08.15.
 */
public class DRQA_All extends BaseExecutable {

    /** The number of RQA compute jobs: #input_files * #embedding_methods * #noise_levels */
    static int results_total;
    /** Report progress ever k seconds. */
    static long report_each_seconds = 10;
    /** Last progress report. */
    static long last_report = System.currentTimeMillis();
    /** Start time of the computation. */
    static long started;

    /** the noise ratios to apply to the original data. */
    public static HashMap<Integer, Double> noiseRatiosById = new HashMap<Integer, Double>();
    static { HashMap<Integer, Double> nr = noiseRatiosById;
//        nr.put(0, 0.0);
//        nr.put(4, 0.1); //
        nr.put(5, 0.2);
//        nr.put(6, 0.4); nr.put(7, 0.8); nr.put(8, 1.6); nr.put(9, 3.2); nr.put(10, 6.4); nr.put(11, 12.8);
    }


    /** Encapsulates the result of an RQA compute job. For simplicity, this is a tab-separated string that can be written directly to disk. */
    private static class Result{
        // the number of produced results so far.
        volatile static int result_id = 1;
        String resultline;
        public Result(String resultline) {
            this.resultline = resultline;
            result_id++;
            if(System.currentTimeMillis()-last_report > report_each_seconds*1000){
                int remaining_jobs = results_total - result_id;
                long average_job_time = (System.currentTimeMillis() - started) / result_id;
                long remaining_time = remaining_jobs * average_job_time / 1000;
                final long elapsed_time_minutes = (last_report - started)/1000/60;
                System.out.println(String.format("finished %.2f%% elapsed time %sh:%sm estimated time left %s minutes %s seconds", 100.0*result_id/results_total, elapsed_time_minutes/60, elapsed_time_minutes, remaining_time / 60, remaining_time%60));
                last_report = System.currentTimeMillis();
            }
        }
    }

    /** the fraction of the maximum phase space trajectory diameter that gives the recurrence threshold. */
    public static double fractionPSTD = 0.05;

    public static void main(String[] args) throws InterruptedException, ExecutionException, IOException {

        if(args.length < 3) {
            System.out.println("Usage: <input_dir> <file_prefix> <out_path> [<num_threads>]");
            System.exit(0);
        }
        String inputDir = args[0];
        String file_prefix = args[1];
        String outPath = args[2];
        int num_threads = args.length > 3 ? Integer.parseInt(args[3]) : Runtime.getRuntime().availableProcessors()-1 ;

        char[] passchars;
        final String jdbc_connection;

        if(outPath.startsWith("jdbc")){
            Console cnsl = System.console();
            passchars = cnsl.readPassword("Database password: ");
            System.out.println("Storing results in database.");
            outPath += "&password=" + new String(passchars);
            jdbc_connection = outPath;
        } else
            jdbc_connection = null;


        System.out.println("Creating Tasks...");

        List<Callable<Result>> tasks = new ArrayList<Callable<Result>>();
        List<File> filteredFiles = getFilteredFiles(inputDir, file_prefix);
        for (int i = 0; i < filteredFiles.size(); i++) {
            final File file = filteredFiles.get(i);
            final double[][] trajectory = BaseExecutable.readFile(file.getAbsolutePath(), ",");
            final String class_name = TimeSeriesGenerator.findClassName(file.getName());
            if(i%100==0) System.out.println(String.format("Reading File %s/%s", i, filteredFiles.size()));

            for (DRQA.EmbeddingMethod embedding_method : new DRQA.EmbeddingMethod[]{DRQA.EmbeddingMethod.ORIGINAL_TRAJECTORY}) {  //DRQA.EmbeddingMethod.ONE_DIMENSION, DRQA.EmbeddingMethod.TIME_DELAY_EMBEDDING_FIRST_MIN_UNBOUNDED, DRQA.EmbeddingMethod.ORIGINAL_TRAJECTORY,
//            for (DRQA.EmbeddingMethod embedding_method : new DRQA.EmbeddingMethod[]{DRQA.EmbeddingMethod.TIME_DELAY_EMBEDDING_FIRST_MIN_UNBOUNDED}) {
//                System.out.println("computing only the TDE TIME_DELAY_EMBEDDING_FIRST_MIN_UNBOUNDED embedding variant!");
//            for (DRQA.EmbeddingMethod embedding_method : DRQA.EmbeddingMethod.values()) {
                final DRQA.EmbeddingMethod embedding_method_task = embedding_method;
                for(Integer noise_ratio_id : noiseRatiosById.keySet()){
                    final int noise_ratio_id_task = noise_ratio_id;
                    Callable<Result> c = new Callable<Result>() {
                        @Override
                        public Result call() throws Exception { return compute(file.getName(), class_name, trajectory, embedding_method_task, noise_ratio_id_task, jdbc_connection); }
                    };
                    tasks.add(c);
                }
            }
        }

        results_total = tasks.size();
        System.out.println(String.format("Processing %s tasks in %s threads.", results_total, num_threads));

        ExecutorService exec = Executors.newFixedThreadPool(num_threads);

        started = System.currentTimeMillis();
        List<Future<Result>> results = exec.invokeAll(tasks);

        // write output only to disk if no database connection was given
        BufferedWriter writer = null;
        if(jdbc_connection == null) {
            writer = Files.newBufferedWriter(Paths.get(outPath), Charset.defaultCharset(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
            writer.write(DRQA.separatedHeader()+"\n");
        }


        try {
            if(jdbc_connection == null) {
                int written_lines = 0;
                for (Future<Result> fr : results) {
                    writer.write(fr.get().resultline);
                    writer.write("\n");
                    written_lines++;
                }
                long elapsed = System.currentTimeMillis() - started;
                System.out.println(String.format("Elapsed time: %sh:%sm", elapsed/1000/60/60, elapsed/1000/60));
                System.out.println(String.format("Lines written: %s", written_lines));
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
         finally {
            exec.shutdown();
            if(jdbc_connection == null) writer.close();
        }
    }

    /**
     * Computes a single RQA job.
     * @param id the id of the original data.
     * @param trajectory the trajectory data.
     * @param embedding_method the offset of the embedding method in the DRQA.EmbeddingMethod.values() array
     * @param noise_ratio_id the offset of the noise ratio (in the noiseRatios array) to be applied to the data (relative to threshold)
     * @return the computed RQA values
     * @throws InterruptedException
     */
    public static Result compute(String id, String class_name, double[][] trajectory, DRQA.EmbeddingMethod embedding_method, int noise_ratio_id, String jdbc_conn) throws InterruptedException {

        double noise_ratio = noiseRatiosById.get(noise_ratio_id);

        double[][] noisy_trajectory = TimeSeriesGenerator.addNoise(trajectory, fractionPSTD, noise_ratio);
        double[][] embedded_trajectory;

        int dimension = trajectory.length;
        int delay = -1;

        // reduce to one dimension, reconstruct if necessary
        switch(embedding_method){
            case ONE_DIMENSION:
                // use only one dimension of the data and do not embed.
                embedded_trajectory = new double[][]{noisy_trajectory[0]};
                dimension = 1;
                break;
            case TIME_DELAY_EMBEDDING:
            case TIME_DELAY_EMBEDDING_SMALLEST_MIN:
            case TIME_DELAY_EMBEDDING_FIRST_MIN_UNBOUNDED:
                PhaseSpaceReconstructed reconstructed = new PhaseSpaceReconstructed(noisy_trajectory[0], 1);
                embedded_trajectory = reconstructed.trajectory;
                dimension = reconstructed.embedding_dimension;
                delay = reconstructed.embedding_delay;
                break;
            case TIME_DELAY_EMBEDDING_DIM2:
                PhaseSpaceReconstructed reconstructed_dim2 = new PhaseSpaceReconstructed(noisy_trajectory[0], 2);
                embedded_trajectory = reconstructed_dim2.trajectory;
                dimension = reconstructed_dim2.embedding_dimension;
                delay = reconstructed_dim2.embedding_delay;
                break;
            default:
                // leave (noisy) data as is
                embedded_trajectory = noisy_trajectory;
                break;
        }

        String result = null;
        Connection connection = null;
        if(jdbc_conn != null){
            try { Class.forName("com.mysql.jdbc.Driver").newInstance(); } catch (Exception ex) { System.err.println("Could not initiate JDBC driver."); ex.printStackTrace(); }
            try { connection = DriverManager.getConnection(jdbc_conn); } catch (SQLException e) { e.printStackTrace(); }
        }

        DRQA drqa = new DRQA(embedded_trajectory, embedded_trajectory, fractionPSTD);
        // with all minimum line lengths set to 2
        drqa.computeRQA(2, 2, 2);
        if(connection != null)
            drqa.writeResultDB(id, class_name, embedding_method.ordinal(), dimension, delay, embedded_trajectory[0].length, DRQA.LMinMethod.L_MIN_FIX.ordinal(), 2, noise_ratio_id, connection);
        else
            result = drqa.toSeparatedString(id, class_name, embedding_method.ordinal(), dimension, delay, DRQA.LMinMethod.L_MIN_FIX.ordinal(), 2, noise_ratio_id, embedded_trajectory[0].length);

        if(connection != null) try { connection.close(); } catch (SQLException e) { e.printStackTrace(); }

        return new Result(result);
    }
}

