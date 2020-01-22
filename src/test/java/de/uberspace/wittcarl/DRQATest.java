package de.uberspace.wittcarl;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import de.uberspace.wittcarl.phasespace.PhaseSpaceDistribution;
import junit.framework.TestCase;

import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.sql.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertArrayEquals;

/**
 * Validates the line length histograms and RQA measures of toy examples.
 */
public class DRQATest extends TestCase {

    // A short example time series
    final double[][] trajectory = new double[][]{ new double[]{0, 0, 2, 3, 5, 5, 5, 3, 2, 0, 0, 2, 3, 0} };

    final double eps = 0.0;

    /**
     * Tests the computation of DRQA line length histograms on a toy example that covers definite and indefinite lines of length 1 and length >1.
     */
    public void testHistogramComputation() {

        DRQA.conditional_ww_limit = 10;
        DRQA drqa = new DRQA(trajectory, eps);

        // compare the computed histograms (diagonal, vertical, ...) against the correct histograms
        assertArrayEquals(new long[]{0, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, drqa.l_hist);
        assertArrayEquals(new long[]{0, 18, 5, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, drqa.v_hist);
        assertArrayEquals(new long[]{0, 0, 8, 3, 3, 3, 0, 5, 0, 0, 0, 0, 0, 0, 0}, drqa.w_hist);
        assertArrayEquals(new long[]{0, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, drqa.r_hist);
        assertArrayEquals(new long[]{0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, drqa.l_hist_indefinite);
        assertArrayEquals(new long[]{0, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, drqa.v_hist_indefinite);
        assertArrayEquals(new long[]{0, 3, 6, 3, 3, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0}, drqa.w_hist_indefinite);
        assertArrayEquals(new long[]{0, 6, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, drqa.r_hist_indefinite);

        double[][] expectedCRT =new double[][]{
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 3, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 3, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 5, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

        assertTrue(Arrays.deepEquals(expectedCRT, drqa.conditional_ww.toArray()));

    }

    /**
     * Tests the computation of line length histogram descriptors.
     */
    public void testHistogramDescriptors() {

        DRQA drqa = new DRQA(trajectory, eps);
        drqa.computeRQA(10,10,10);  // should not change the result.
        drqa.computeRQA(2, 2, 2);

        Map<DRQA.HistogramStatistic, Double>
                expected_diagonal_rqa = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_vertical_rqa = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_white_vertical_rqa = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_orthogonal_rqa = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_diagonal_rqa_definite = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_vertical_rqa_definite = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_white_vertical_rqa_definite = new HashMap<DRQA.HistogramStatistic, Double>(),
                expected_orthogonal_rqa_definite = new HashMap<DRQA.HistogramStatistic, Double>();

        expected_diagonal_rqa.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 0.6365141682948128);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.SORTEDNESS, 0.6666666666666666);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.MAX, 14.0);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.DEFINITENESS, 0.2727272727272727);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.ENTROPY, 1.0986122886681096);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.FILTER_RATIO, 0.6060606060606061);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.SPARSENESS, 0.8333333333333334);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.07142857142857142);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 3.0/7.0);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.AVERAGE, 6.666666666666666);
        expected_diagonal_rqa.put(DRQA.HistogramStatistic.MEDIAN, 4.);

        expected_vertical_rqa.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 0.0);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.SORTEDNESS, 1.0);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.MAX, 3.0);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.DEFINITENESS, 0.7115384615384616);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.ENTROPY, 0.5402041423888608);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.FILTER_RATIO, 0.5576923076923077);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.SPARSENESS, 0.0);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.3333333333333333);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 1.0);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.AVERAGE, 2.230769230769231);
        expected_vertical_rqa.put(DRQA.HistogramStatistic.MEDIAN, 2.);

        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 1.3296613488547582);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.SORTEDNESS, 0.65);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.MAX, 7.0);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.DEFINITENESS, 0.6041666666666666);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.ENTROPY, 1.6470013963439956);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.FILTER_RATIO, 1.0);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.SPARSENESS, 0.16666666666666666);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.14285714285714285);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 2.0/4.0);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.AVERAGE, 3.5999999999999996);
        expected_white_vertical_rqa.put(DRQA.HistogramStatistic.MEDIAN, 3.);

        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 0.0);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.SORTEDNESS, 0.5);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.MAX, 5.0);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.DEFINITENESS, 0.42105263157894735);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.ENTROPY, 0.6931471805599453);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.FILTER_RATIO, 0.42105263157894735);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.SPARSENESS, 0.6666666666666666);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.2);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 1.0);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.AVERAGE, 4.0);
        expected_orthogonal_rqa.put(DRQA.HistogramStatistic.MEDIAN, 3.);

        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 0.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.SORTEDNESS, 1.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.MAX, 2.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.DEFINITENESS, 0.2727272727272727);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.ENTROPY, 0.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.FILTER_RATIO, 0.6060606060606061);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.SPARSENESS, 0.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.5);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 1.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.AVERAGE, 2.0);
        expected_diagonal_rqa_definite.put(DRQA.HistogramStatistic.MEDIAN, 2.);

        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 0.0);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.SORTEDNESS, 1.0);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.MAX, 3.0);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.DEFINITENESS, 0.7115384615384616);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.ENTROPY, 0.6615632381579821);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.FILTER_RATIO, 0.5576923076923077);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.SPARSENESS, 0.0);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.3333333333333333);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 1.0);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.AVERAGE, 2.375);
        expected_vertical_rqa_definite.put(DRQA.HistogramStatistic.MEDIAN, 2.);

        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 1.3296613488547582);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.SORTEDNESS, 0.6363636363636364);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.MAX, 7.0);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.DEFINITENESS, 0.6041666666666666);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.ENTROPY, 1.5196682491027627);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.FILTER_RATIO, 1.0);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.SPARSENESS, 0.3333333333333333);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.14285714285714285);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 0.5);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.AVERAGE, 3.954545454545454);
        expected_white_vertical_rqa_definite.put(DRQA.HistogramStatistic.MEDIAN, 3.);

        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.DIFFERENCE_ENTROPY, 0.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.SORTEDNESS, 0.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.MAX, 3.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.DEFINITENESS, 0.42105263157894735);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.ENTROPY, 0.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.FILTER_RATIO, 0.42105263157894735);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.SPARSENESS, 1.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.INVERSE_MAX, 0.3333333333333333);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.LOCAL_MAXIMA, 1.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.AVERAGE, 3.0);
        expected_orthogonal_rqa_definite.put(DRQA.HistogramStatistic.MEDIAN, 3.);

        assertEquals(expected_diagonal_rqa, drqa.diagonal_rqa);
        assertEquals(expected_vertical_rqa, drqa.vertical_rqa);
        assertEquals(expected_white_vertical_rqa, drqa.white_vertical_rqa);
        assertEquals(expected_orthogonal_rqa, drqa.orthogonal_rqa);
        assertEquals(expected_diagonal_rqa_definite, drqa.diagonal_rqa_definite);
        assertEquals(expected_vertical_rqa_definite, drqa.vertical_rqa_definite);
        assertEquals(expected_white_vertical_rqa_definite, drqa.white_vertical_rqa_definite);
        assertEquals(expected_orthogonal_rqa_definite, drqa.orthogonal_rqa_definite);

        assertEquals(Math.max(2, Math.round(37.0 / 26.)), drqa.estimate_l_min(), 1e-10);

    }

    public void testWriteDB() throws SQLException {
        try { Class.forName("com.mysql.jdbc.Driver").newInstance(); } catch (Exception ex) { System.err.println("Could not initiate JDBC driver."); ex.printStackTrace(); }

        try {
            Connection conn = DriverManager.getConnection("jdbc:mysql://localhost:3306/wittcarl_recurrence_plot_clustering?user=root&password=mysql");
            TimeSeriesGenerator tsg = new TimeSeriesGenerator(1,1000, "lorenz");
            double[][] trajectory = tsg.lorenz_standard[0];
            DRQA drqa = new DRQA(trajectory, trajectory, 0.05);
            drqa.computeRQA(2, 2, 2);
            drqa.writeResultDB("IN_MEMORY", "lorenz", DRQA.EmbeddingMethod.ORIGINAL_TRAJECTORY.ordinal(), 3, -1, trajectory[0].length, DRQA.LMinMethod.L_MIN_FIX.ordinal(), 2, 0, conn);

        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
        }
    }

    public void testSerialize() throws IOException, ClassNotFoundException {

        SparseDoubleMatrix2D matrix2D = new SparseDoubleMatrix2D(new double[][]{new double[]{1,2,3}, new double[]{4,5,6}, new double[]{7,8,9}});
        System.out.println(String.format("matrix2D: %s", matrix2D));
        try { Class.forName("com.mysql.jdbc.Driver").newInstance(); } catch (Exception ex) { System.err.println("Could not initiate JDBC driver."); ex.printStackTrace(); }

        try {
            Connection conn = DriverManager.getConnection("jdbc:mysql://localhost:3306/wittcarl_recurrence_plot_clustering?user=root&password=mysql");
            PreparedStatement preps = conn.prepareStatement("INSERT INTO `wittcarl_recurrence_plot_clustering`.`blob` (`id`, `signature`) VALUES (null, ?);");
            preps.setObject(1, matrix2D);
            preps.execute();

            PreparedStatement get = conn.prepareStatement("SELECT id, signature FROM `blob` WHERE id=3");
            ResultSet resultSet = get.executeQuery();
            resultSet.next();
            int id = resultSet.getInt(1);

            InputStream is = resultSet.getBlob(2).getBinaryStream();
            ObjectInputStream oip = new ObjectInputStream(is);
            Object object = oip.readObject();
            String className = object.getClass().getName();
            oip.close();
            is.close();
            resultSet.close();

            // de-serialize list a java object from a given objectID
            SparseDoubleMatrix2D restored = (SparseDoubleMatrix2D) object;

            System.out.println(String.format("id: %s", id));
            System.out.println(String.format("restored: %s", restored));
            conn.close();

        } catch (SQLException ex) {
            // handle any errors
            System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
        }
    }

    public void testCompress(){

        SparseDoubleMatrix2D matrix = new SparseDoubleMatrix2D(new double[][]{
                new double[]{0,0,0,0},            // 0 0 0 0
                new double[]{0,1,10,100},         // 0 0 1 2
                new double[]{0,100000,100,1}      // 0 3 2 0
        });

        for(short[] row : DRQA.compress_hist(matrix)){
            System.out.println(Arrays.toString(row));
        }
    }

    /**
     * Compare approximate and exact diameter computation.
     * Approximate is very fast but provides only a rough upper bound.
     * For a discretized pairwise distance distribution, this is absolutely sufficient.
     */
    public void testMaxPhaseSpaceDiameter(){

        int tsLen = 25000;
        TimeSeriesGenerator tsg = new TimeSeriesGenerator(1, tsLen, "lorenz");
        double[][][][] trajectories = tsg.getAllTrajectories();
        double[][] ts = trajectories[4][0];
        System.out.println("Computing max phase space diameter.");
        long tic = System.currentTimeMillis();
        double diam = PhaseSpaceDistribution.maxPhaseSpaceDiameter(ts, ts);
        System.out.println("diam = " + diam);
        long toc = System.currentTimeMillis();
        System.out.println(toc-tic);

        System.out.println("Computing max phase space diameter, approximate.");
        tic = System.currentTimeMillis();
        double diamApprox = PhaseSpaceDistribution.maxPhaseSpaceDiameterApproximate(ts);
        System.out.println("diamApprox = " + diamApprox);
        toc = System.currentTimeMillis();
        System.out.println(toc-tic);
        System.out.println("Relative error: " + (Math.abs(diamApprox - diam)/diam));

//        GeometryFactory factory = new GeometryFactory(new PrecisionModel(PrecisionModel.FIXED));
//        Point[] points = new Point[tsLen];
//        for (int i = 0; i < tsLen; i++) {
//            Coordinate[] coordinates = new Coordinate[ts.length];
//            for (int j = 0; j < ts.length; j++) {
//                coordinates[j] = new Co
//            }
//            points[i] = new Point(new CoordinateArraySequence(coordinates), factory);
//        }
//        MultiPoint data = new MultiPoint(points, factory);
//        MinimumDiameter minimumDiameter = new MinimumDiameter(geom);

    }

    public void testPhaseSpaceApproxDistribution(){

        int tsLen = 500;
        TimeSeriesGenerator tsg = new TimeSeriesGenerator(1, tsLen, null);
        double[][][][] trajectories = tsg.getAllTrajectories();
        for (int systemIdx = 0; systemIdx < 9; systemIdx++) {
            double[][] ts = trajectories[systemIdx][0];
            long[] dist = PhaseSpaceDistribution.approximateDistanceDistribution(ts);
            System.out.println(Arrays.toString(dist));
        }

    }

}