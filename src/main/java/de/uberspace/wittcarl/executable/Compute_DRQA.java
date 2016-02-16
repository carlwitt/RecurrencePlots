package de.uberspace.wittcarl.executable;

import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;

import java.io.File;
import java.nio.file.*;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 12.08.15.
 */
public class Compute_DRQA extends BaseExecutable {



    public static void main(String[] args) throws Exception {

        if (args.length < 4) {
            System.err.println("Usage: input_dir file_prefix recurrence_threshold_fraction_max_phase_space_diameter out_file [crt_limit] [write_rp]");
            System.exit(0);
        }

        String path = args[0];
        String prefix = args[1];
        double fracMaxPSDDiameter = Double.parseDouble(args[2]);
        String outfile = args[3];
        int crtLimit = args.length > 4 ? Integer.parseInt(args[4]) : 0;
        boolean writeRP = args.length > 5 && Boolean.parseBoolean(args[5]);


        DRQA.conditional_ww_limit = crtLimit;

        System.out.println(String.format("Input Directory: %s Prefix Filter: %s, data points\nThreshold: %s", path, prefix, fracMaxPSDDiameter));

        String header = DRQA.separatedHeader();
        if(! outfile.equals("stdout"))
            Files.write(Paths.get(outfile), (header+"\n").getBytes(), StandardOpenOption.CREATE);

        for(File file : getFilteredFiles(path, prefix))
            processFile(file.getAbsolutePath(), outfile, fracMaxPSDDiameter, writeRP, crtLimit);

    }

    private static void processFile(String filename, String outputfile, double fracMaxPSDDiameter, boolean writeRP, int crtLimit) throws Exception {

        System.out.println(String.format("Processing: %s", filename));
        double[][] values = BaseExecutable.readFile(filename, ",");

        // record start time
        long started = System.currentTimeMillis();

        DRQA drqa = new DRQA(values, values, fracMaxPSDDiameter);
        System.out.println(String.format("DRQA computed in: %s ms", System.currentTimeMillis() - started));
        drqa.computeRQA(2, 2, 2);

        if(outputfile.equals("stdout"))
            drqa.print();
        else{
            String line = drqa.toSeparatedString(filename, TimeSeriesGenerator.findClassName(filename), DRQA.EmbeddingMethod.ORIGINAL_TRAJECTORY.ordinal(), values.length, -1, DRQA.LMinMethod.L_MIN_FIX.ordinal(), 2, 0, values[0].length) + "\n";
            Files.write(Paths.get(outputfile), line.getBytes(), StandardOpenOption.APPEND);
        }

        Files.write(Paths.get(filename+"_CRT.tsv"), drqa.conditionalRecurrenceTsv(drqa.conditional_ww, crtLimit).toString().getBytes());
        drqa.writeCRTImages(filename);

        if(writeRP) drqa.writeRP(values, values, fracMaxPSDDiameter, filename+"_RP.png");

    }

}
