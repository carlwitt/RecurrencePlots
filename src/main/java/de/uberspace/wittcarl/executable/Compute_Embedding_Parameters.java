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

import de.uberspace.wittcarl.phasespace.PhaseSpaceReconstructed;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

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

@Deprecated
public class Compute_Embedding_Parameters extends BaseExecutable {

    public static void main(String[] args) throws IOException {

        if (args.length < 3) {
            System.err.println("Usage: input_dir file_prefix rqa_results_file [num_processes]");
            System.exit(0);
        }

        String input_dir = args[0];
        String prefix = args[1];
        final String outputfile = args[2];
        int num_threads = args.length > 3 ? Integer.parseInt(args[3]) : Runtime.getRuntime().availableProcessors()-1 ;

        // integer given, interpret as dimension
        System.out.println(String.format("Input Directory: %s", input_dir));
        System.out.println(String.format("File Prefix: %s", prefix));

        System.out.println("Creating Tasks...");

        List<Callable<String>> tasks = new ArrayList<Callable<String>>();
        List<File> filteredFiles = getFilteredFiles(input_dir, prefix);
        for (final File file : filteredFiles) {
            Callable<String> c = new Callable<String>() { @Override public String call() {
                try { return processFile(file); }
                catch (IOException e) { e.printStackTrace(); return null; }
            } };
            tasks.add(c);
        }

        ExecutorService exec = Executors.newFixedThreadPool(num_threads);
        System.out.println(String.format("Processing %s tasks in %s threads.", tasks.size(), num_threads));

        BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputfile), Charset.defaultCharset(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
        writer.write("path\tclass\tdim\tsmallest_local_min\tfirst_local_min\n");

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
            System.out.println(String.format("Lines written: %s", written_lines));
        } catch (Exception e) {
            e.printStackTrace();
        }
        finally {
            exec.shutdown();
            writer.close();
        }


    }

    private static String processFile(File file) throws IOException {

        final String absolutePath = file.getAbsolutePath();
        final String filename = file.getName();

        System.out.println(String.format("Processing: %s", absolutePath));

        double[][] trajectory = BaseExecutable.readFile(absolutePath, ",");
        double[][] noisyTrajectory = TimeSeriesGenerator.addNoise(trajectory, 0.05, 0.1);

        PhaseSpaceReconstructed phaseSpaceReconstructed = new PhaseSpaceReconstructed(noisyTrajectory[0], 1);

        return String.format("%s\t%s\t%s\t%s\t%s",filename, filename.split("\\.")[0].subSequence(6,filename.split("\\.")[0].length()-2), phaseSpaceReconstructed.embedding_dimension, phaseSpaceReconstructed.embedding_delay, phaseSpaceReconstructed.min_embedding_delay);
    }

}
