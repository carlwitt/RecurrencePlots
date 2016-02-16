package de.uberspace.wittcarl.executable;

import junit.framework.TestCase;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 20.08.15.
 */
public class DRQAConvergenceTest extends TestCase {

    public void testComputeConvergenceTimes() throws Exception {

        DRQAConvergence.computeConvergenceTimes(100, 140, "test-convergence.tsv", 0.01, 100);

    }
}
