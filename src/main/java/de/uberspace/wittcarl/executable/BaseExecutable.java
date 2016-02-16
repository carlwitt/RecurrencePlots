package de.uberspace.wittcarl.executable;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.MalformedInputException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 17.08.15.
 *
 * Provides basic in/out methods and file system access.
 *
 */
public class BaseExecutable {

    public static List<File> getFilteredFiles(String path, String prefix) {
        List<File> filenames = new ArrayList<File>();
        final File[] files = new File(path).listFiles();
        if (files != null)
            for (File file : files)
                if (!file.isDirectory() && file.getName().startsWith(prefix) && !file.getName().startsWith(".")) filenames.add(file);
        return filenames;
    }

    public static void writeStringToFile(String path, String content) {
        try {
            Files.write(Paths.get(path), content.getBytes());
        } catch (IOException e) {
            System.err.println("Couldn't write file " + path);
            e.printStackTrace();
        }
    }

    public static double[][] readFile(String path, String delimiter) {

        try {

            List<String> lines;
            Path file = Paths.get(path);
            try {
                lines = Files.readAllLines(file, Charset.defaultCharset());
            } catch (MalformedInputException e) {
                System.err.println("File is not a text file, skipping: " + path);
                return new double[][]{};
            }
            // remove trailing empty lines
            for (int i = lines.size()-1; i >= 0; i--) {
                if(lines.get(i).length() == 0) lines.remove(i);
            }

            int trajectory_dimensionality = lines.get(0).split(",").length;
            double[][] values = new double[trajectory_dimensionality][lines.size()];

            for (int i = 0; i < lines.size(); i++) {
                String row = lines.get(i);
                final String[] lineValues = row.split(delimiter);
                for (int j = 0; j < trajectory_dimensionality; j++)
                    values[j][i] = Double.parseDouble(lineValues[j]);
            }

            return values;

        } catch (IOException e) {
            System.err.println("Could not read file " + path);
            return new double[][]{};
        }

    }

}
