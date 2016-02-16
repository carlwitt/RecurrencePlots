package de.uberspace.wittcarl.executable;

import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import junit.framework.TestCase;

import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 17.08.15.
 */
public class BaseExecutableTest extends TestCase{

    public void testReadFile() throws Exception {

        final String file = "data/junit-test-data/lorenz-noisy-10000/lorenz-10000-0.0.csv";
        double[][] trajectory = BaseExecutable.readFile(file, ",");
        assertEquals(TimeSeriesGenerator.trajectoryToCsvString(trajectory).toString(), new String(Files.readAllBytes(Paths.get(file))));

    }

}