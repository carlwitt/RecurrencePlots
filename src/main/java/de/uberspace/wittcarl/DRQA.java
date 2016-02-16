package de.uberspace.wittcarl;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.jet.random.Empirical;
import cern.jet.random.engine.DRand;
import de.uberspace.wittcarl.phasespace.PhaseSpaceDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * Created by Carl Witt on 01.08.15.
 * Computes the line length histograms of RPs of space trajectories, writes RPs and CRT plots.
 * Uses a database or textfiles to output the results.
 * Provides bootstrapping facilities for line length histograms.
 */
public class DRQA {

    public static boolean CRT_LOG_SCALE = false;
    long recurrence_points;
    double recurrence_rate;
    public final double recurrence_threshold;

    public long[] l_hist;               // diagonal lines
    public long[] v_hist;               // vertical lines
    public long[] w_hist;               // white vertical lines
    public long[] r_hist;               // reverse diagonal lines
    public long[] l_hist_indefinite;    // diagonal lines that touch the borders of the matrix, we have only a lower bound on their length
    public long[] r_hist_indefinite;    // one could count simply their number, but for determinism and laminarity, we can use the lower bound
    public long[] v_hist_indefinite;    // to incorporate the pixels on indefinite lines according to the lower bound of their length
    public long[] w_hist_indefinite;
    // distribution of pairs of subsequent white vertical lines
    static public int conditional_ww_limit = 1000;
    public SparseDoubleMatrix2D conditional_ww;

    public String crtStatistics() {

        double[][] histo = conditional_ww.toArray();
        double totalWeight = conditional_ww.zSum();

        // mean
        double mean_row = 0, mean_col = 0;
        for (int row = 0; row < histo.length; row++) {
            for (int col = 0; col < histo[row].length; col++) {
                mean_row += histo[row][col] * row;
                mean_col += histo[row][col] * col;
            }
        }
        mean_row /= totalWeight;
        mean_col /= totalWeight;

        // variance
        double var_row = 0, var_col = 0;
        for (int row = 0; row < histo.length; row++) {
            for (int col = 0; col < histo[row].length; col++) {
                var_row += histo[row][col] * (row - mean_row) * (row - mean_row);
                var_col += histo[row][col] * (col - mean_col) * (col - mean_col);
            }
        }
        var_row /= totalWeight;
        var_col /= totalWeight;

        // covar row col
        double correlation = 0;
        for (int row = 0; row < histo.length; row++) {
            for (int col = 0; col < histo[row].length; col++) {
                correlation += histo[row][col] * (row-mean_row) * (col-mean_col);
            }
        }
        correlation /= totalWeight * Math.sqrt(var_row) * Math.sqrt(var_col);

        // global maximum ...
        int max_row = 0, max_col = 0;
        double max_val = histo[0][0];
        for (int row = 3; row < histo.length; row++) {
            for (int col = 3; col < histo[row].length; col++) {
                if(histo[row][col] > max_val){
                    max_row = row;
                    max_col = col;
                    max_val = histo[row][col];
                }
            }
        }
        // ... and number of local maxima
        int local_maxima = 0;
        int neighborhoodsize = 5;
        for (int row = neighborhoodsize; row < histo.length-neighborhoodsize; row++) {
            for (int col = neighborhoodsize; col < histo[row].length-neighborhoodsize; col++) {

                double current_value = histo[row][col];

                // check neighborhood
                boolean largerThanAll = true;
                for (int row_offset = -neighborhoodsize; (row_offset <= neighborhoodsize) && largerThanAll; row_offset++) {
                    for (int col_offset = -neighborhoodsize; (col_offset <= neighborhoodsize) && largerThanAll; col_offset++) {
                        if(row_offset == 0 && col_offset == 0) continue;    // do not compare to self
                        if(histo[row+row_offset][col+col_offset] > current_value) largerThanAll = false;
                    }
                }
                if(largerThanAll) local_maxima++;
            }
        }

        // entropy
        double entropy = 0;
        if(totalWeight > 0) {
            for (double[] row : histo) {
                for (double bin : row) {
                    if(bin < 1) continue;
                    final double prob = bin / totalWeight;
                    entropy += prob * Math.log(prob);
                }
            }
        } else { entropy = 0; }
        entropy *= -1;

        return String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",mean_row, mean_col, correlation, max_row, max_col, local_maxima, entropy);
    }

    public String printableString(HistogramStatistic[] measure_selection) {
        double[] allMeasures = allRQAMeasures(measure_selection);
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < allMeasures.length; i++) {
            HistogramStatistic h = measure_selection[i%measure_selection.length];
            b.append(LineType.values()[i/measure_selection.length].name());
            b.append(h.name());
            b.append(": ");
            b.append(String.format("%.2f", allMeasures[i]));
            b.append("\n");
        }
        return b.toString();
    }
//    static public int conditional_vv_limit = 1000;
//    public SparseDoubleMatrix2D conditional_vv;
//    static public int conditional_wv_limit = 1000;
//    public SparseDoubleMatrix2D conditional_wv;

    public enum LineType{ DIAGONAL, VERTICAL, WHITE_VERTICAL, ORTHOGONAL }
    public enum HistogramStatisticKind{ RQA, DRQA }
    public enum HistogramStatistic{ FILTER_RATIO, AVERAGE, MEDIAN, MAX, INVERSE_MAX, ENTROPY, DEFINITENESS, SORTEDNESS, LOCAL_MAXIMA, SPARSENESS, DIFFERENCE_ENTROPY, }
    public static final HistogramStatistic[] STANDARD_RQA = new HistogramStatistic[]{HistogramStatistic.FILTER_RATIO, HistogramStatistic.AVERAGE, HistogramStatistic.MAX, HistogramStatistic.INVERSE_MAX, HistogramStatistic.ENTROPY};
    public enum CRTStatistic { MEAN_ROW, MEAN_COL, CORRELATION, MAX_ROW, MAX_COL, LOCAL_MAXIMA, ENTROPY }

    public HashMap<HistogramStatistic, Double>
            diagonal_rqa = new HashMap<HistogramStatistic, Double>(),
            vertical_rqa = new HashMap<HistogramStatistic, Double>(),
            white_vertical_rqa = new HashMap<HistogramStatistic, Double>(),
            orthogonal_rqa = new HashMap<HistogramStatistic, Double>(),
            diagonal_rqa_definite = new HashMap<HistogramStatistic, Double>(),
            vertical_rqa_definite = new HashMap<HistogramStatistic, Double>(),
            white_vertical_rqa_definite = new HashMap<HistogramStatistic, Double>(),
            orthogonal_rqa_definite = new HashMap<HistogramStatistic, Double>();
    HashMap<CRTStatistic, Double> crt_statistics = new HashMap<CRTStatistic, Double>();

    /**
     * @param trajectory the phase space trajectory, where trajectory[j][i] refers to the j-th component (dimension) of the i-th phae space trajectory point
     * @param eps the threshold.
     */
    public DRQA(double[][] trajectory, double eps){

        recurrence_threshold = eps;

        computeHistograms(trajectory, trajectory, eps);
        computeRQA(2, 2, 2);

        // test that all recurrence points sum up nicely
        checkHistogramSums(trajectory[0].length);

    }

    /**
     * @param embedded1 the first phase space trajectory, where embedded[j][i] refers to the j-th component (dimension) of the i-th phae space trajectory point
     * @param embedded2 as embedded1 but for cross recurrence plots.
     * @param fractionMaxPSTDiameter The fraction of the maximum phase space diameter to use as recurrence threshold. A common value in the literature is 0.05.
     */
    public DRQA(double[][] embedded1, double[][] embedded2, double fractionMaxPSTDiameter){

        recurrence_threshold = fractionMaxPSTDiameter * PhaseSpaceDistribution.maxPhaseSpaceDiameter(embedded1, embedded2);

        computeHistograms(embedded1, embedded2, recurrence_threshold);

        // test that all recurrence points sum up nicely
        checkHistogramSums(embedded1[0].length);

    }

    public static long[] approximateDistanceDistribution(double[][] ts){
        return approximateDistanceDistribution(ts, 500, 1.0);
    }

    /**
     * Approximates the distribution of pairwise distances between phase space points. The approximation consists in quantizing distances.
     * Relative to the approximate maximum phase space diameter {@link #maxPhaseSpaceDiameterApproximate(double[][])}
     * @param ts The trajectory [dimension][time step]
     * @param numBins The number of different distances.
     * @param upTo The maximum distance to record, relative to the approx. max. phase space diameter.
     * @return
     */
    public static long[] approximateDistanceDistribution(double[][] ts, int numBins, double upTo){
        double approxMax = upTo * DRQA.maxPhaseSpaceDiameterApproximate(ts);
        final double scale = numBins / approxMax;
        long[] distanceCount = new long[numBins+1];
        final int dim = ts.length;
        final int N = ts[0].length;

        for (int col = 0; col < N; col++) {
            for (int row = 0; row < N; row++) {
                // compute euclidean distance in dim dimensions
                double squared_distance = 0;
                for (int d = 0; d < dim; d++) {
                    squared_distance += Math.pow(ts[d][col] - ts[d][row], 2);// * (ts[d][col] - ts[d][row]);
                }
                int binIdx = (int) (scale * Math.sqrt(squared_distance));
                if(binIdx <= numBins) distanceCount[binIdx]++;
            }
        }
        distanceCount[0] = Double.doubleToLongBits(approxMax);
        return distanceCount;

    }

    /**
     * A linear space algorithm to compute the line length histograms and empirical conditional recurrence time distribution.
     * @param embedded1 The first phase space trajectory where embedded[i][j] refers to the i-th dimension of the j-th phase space trajectory point.
     * @param embedded2 The second phase space trajectory. As embedded1. For a normal RP, pass the same trajectory as embedded1 and embedded2.
     * @param eps The recurrence threshold (distances are computed according to Euclidean Norm).
     */
    private void computeHistograms(double[][] embedded1, double[][] embedded2, double eps) {// the phase space dimensionality

        final int dim = embedded1.length;
        assert embedded1.length == embedded2.length : "Both trajectories must have the same dimensionality.";

        // sqrt(distance^2) < eps <=> distance^2 < eps^2
        // comparing squared distances to the squared threshold is faster than computing the root of all distances just to get the equivalent result
        double squared_threshold = eps * eps;

        // rqa measures
        recurrence_points = 0;

        final int N = embedded1[0].length;
        assert embedded1[0].length == embedded2[0].length : "Both trajectories must have the same number of observations.";

        l_hist = new long[N+1];
        r_hist = new long[N+1];
        v_hist = new long[N+1];
        w_hist = new long[N+1];
        l_hist_indefinite = new long[N+1];
        r_hist_indefinite = new long[N+1];
        v_hist_indefinite = new long[N+1];
        w_hist_indefinite = new long[N+1];

        // the conditional recurrence P(w|x) denotes the empirical probability to find a white vertical line of length w when the previous white vertical line had length x.
        if(conditional_ww_limit > 0)
            conditional_ww = new SparseDoubleMatrix2D(conditional_ww_limit, conditional_ww_limit);
//        conditional_vv = new SparseDoubleMatrix2D(conditional_vv_limit, conditional_vv_limit);
//        conditional_wv = new SparseDoubleMatrix2D(conditional_wv_limit, conditional_wv_limit);

        // -------------------------------------------------------------------------------------------------------------
        // compute line length histograms
        // -------------------------------------------------------------------------------------------------------------

        // l_carry_previous[i] refers to the line length adjacent to the pixel in row i of the previous column
        // l_carry[i] refers to the line length adjacent to the pixel in row i of the current column
        // having just one is not enough because when updating one overwrites the value (i) that is needed in the next iteration.
        // the carry variables are declared int because the maximum line length is the length of the time series. Due to quadratic runtime time series can't be too long and this won't overflow.
        int[] l_carry_previous = new int[N];
        int[] l_carry = new int[N];
        int[] r_carry_previous = new int[N+1];
        int[] r_carry = new int[N+1];

        int[] temp_carry; // for swapping old and new carry.

        double squared_distance;                                   // sum of squares of vector components, used in the euclidean distance calculation

        // compute the RP in column-major order. requires linear space since the matrix is never stored.
        for (int col = 0; col < N; col++) {

            int last_w = -1;                // the length of the previous white vertical line in this column. -1 denotes that there was none yet.
//            int last_v = -1;                // the length of the previous vertical line in this column. -1 denotes that there was none yet.
            int v_carry = 0;                // the length of the current vertical line
            int w_carry = 0;                // the length of the current white vertical line

            // if a reverse diagonal touches the matrix border, register it with indefinite length
            if(r_carry_previous[0] > 0) r_hist_indefinite[r_carry_previous[0]]++;

            for (int row = 0; row < N; row++) {

                // compute euclidean distance in dim dimensions
                squared_distance = 0;
                for (int d = 0; d < dim; d++) squared_distance += (embedded1[d][col] - embedded2[d][row]) * (embedded1[d][col] - embedded2[d][row]);

                // branch on whether we have a recurrence point or not
                if(squared_distance <= squared_threshold){ // recurrence point

                    recurrence_points++;

                    // diagonal lines (grow)
                    if(row > 0) l_carry[row] = l_carry_previous[row-1] + 1;          // extends the line that ends at (row-1, col-1).
                    else l_carry[0] = 1;

                    // vertical lines (grow)
                    v_carry++;

                    // white vertical lines (finish)
                    if(w_carry > 0){
                        final boolean w_start_known = w_carry < row;   // if the length of the current white vertical line equals the number of previously processed cells, it touches the matrix border and thus it's length is indefinite
                        if(w_start_known){
                            w_hist[w_carry]++;
                            if(last_w > 0 && conditional_ww_limit > 0) increase_sparse_matrix(conditional_ww, last_w, w_carry, 1);
                            last_w = w_carry;
                        } else {
                            w_hist_indefinite[w_carry]++;
                        }
                        w_carry = 0;
                    }

                    // reverse diagonal lines (grow)
                    // the referenced cell must be above the main diagonal (row < col), to continue the line
                    if(row+1 < col-1) r_carry[row] = r_carry_previous[row+1]+1;
                        // if the referenced cell is on the main diagonal, start a new line
                    else if (row+1-(col-1) <= 1) r_carry[row] = 1;
                    // if the referenced cell is below the diagonal, do nothing


                } else { // no recurrence point

                    // diagonal lines (finish)
                    if(row > 0 && l_carry_previous[row-1] > 0){
                        int line_length = l_carry_previous[row-1];
                        int start_row = row - line_length;      // start row is the zero based index of the trajectory point involved in the recurrence point that marks the start the diagonal line
                        int start_col = col - line_length;      // start col is the zero based index of the other trajectory point involved in the recurrence point that marks the start the diagonal line
                        // check whether the line begins at the matrix border. it doesn't end at a matrix border, because we're looking at the white point that ends it right now.
                        final boolean l_start_known = start_row > 0 && start_col > 0;
                        // only above main diagonal
                        if(row < col) if(l_start_known) l_hist[line_length]++; else l_hist_indefinite[line_length]++;
                        l_carry[row] = 0;
                    }

                    // vertical lines (finish)
                    if(v_carry > 0){
                        final boolean v_start_known = v_carry < row;   // if the length of the current vertical line equals the number of readily processed cells, it touches the matrix border and thus it's length is indefinite
                        if(v_start_known){
                            v_hist[v_carry]++;
//                            if(last_w > 0) increase_sparse_matrix(conditional_wv, last_w, v_carry, 1);
//                            if(last_v > 0) increase_sparse_matrix(conditional_vv, last_v, v_carry, 1);
//                            last_v = v_carry;
                        } else v_hist_indefinite[v_carry]++;
                        v_carry = 0;
                    }

                    // white vertical lines (grow)
                    w_carry++;

                    // reverse diagonal lines (finish)
                    if(r_carry_previous[row+1] > 0) {
                        r_hist[r_carry_previous[row+1]]++;  // since we have here the white cell that ends the line (the other side is always ended by the main diagonal), this length is always definite
                        r_carry[row] = 0;
                    }

                } // end branch on recurrence point or not


            } // end for rows (end of current column)

            // the current diagonal line carry will be the previous for the next iteration
            temp_carry = l_carry_previous;
            l_carry_previous = l_carry;
            l_carry = temp_carry; // this is just to reuse the memory since creating and collecting arrays is expensive.
            Arrays.fill(l_carry, 0);
            // the same for the reverse diagonal lines
            temp_carry = r_carry_previous;
            r_carry_previous = r_carry;
            r_carry = temp_carry;
            Arrays.fill(r_carry, 0);

            // one could register here the lines that are finished by means of the end of the matrix, but this wouldn't be correct, since one doesn't know whether they continue or not.
            // instead, we count the lines without registering their length
            if(v_carry > 0) v_hist_indefinite[v_carry]++;
            if(w_carry > 0) w_hist_indefinite[w_carry]++;

        } // end of all columns

        // register the diagonal lines that were ended by the matrix border as lines of indefinite length.
        for (int row = 1; row < N+1; row++) {
            if(l_carry_previous[row-1] > 0){
                int line_length = l_carry_previous[row-1];
                // all these lines are of indefinite length since we are looking at the right margin of the matrix.
                l_hist_indefinite[line_length]++;
            }
        }
        // register all the reverse diagonal lines that are ended by the matrix border (thus, have indefinite length)
        for (int i = 0; i < N - 1; i++) if(r_carry_previous[i] > 0) r_hist_indefinite[r_carry_previous[i]]++;

        recurrence_rate = (double) recurrence_points / N / N;
    }

    private static void increase_sparse_matrix(SparseDoubleMatrix2D matrix2D, int row, int col, double value) {
        double previous = matrix2D.getQuick(row, col);
        matrix2D.setQuick(row, col, previous + value);
    }

    public void writeRP(double[][] trajectory1, double[][] trajectory2, double fractionMaxPSTDiameter, String outpath) {

        BufferedImage off_Image = getRPImage(trajectory1, trajectory2, fractionMaxPSTDiameter);
        File outputfile = new File(outpath);
        try { ImageIO.write(off_Image, "png", outputfile); } catch (IOException e) { e.printStackTrace(); }

    }

    public static BufferedImage getRPImage(double[][] trajectory1, double[][] trajectory2, double fractionMaxPSTDiameter) {
        double eps = fractionMaxPSTDiameter * PhaseSpaceDistribution.maxPhaseSpaceDiameter(trajectory1, trajectory2);

        // the phase space dimensionality
        final int dim = trajectory1.length;
        assert trajectory1.length == trajectory2.length : "Both trajectories must have the same dimensionality.";

        // sqrt(distance^2) < eps <=> distance^2 < eps^2
        // comparing squared distances to the squared threshold is faster than computing the root of all distances just to get the equivalent result
        double squared_threshold = eps * eps;

        final int N = trajectory1[0].length;
        assert trajectory1[0].length == trajectory2[0].length : "Both trajectories must have the same number of observations.";

        double squared_distance;                                   // sum of squares of vector components, used in the euclidean distance calculation

        BufferedImage off_Image = new BufferedImage(N, N, BufferedImage.TYPE_INT_ARGB);

        // compute the RP in column-major order. requires linear space since the matrix is never stored.
        for (int col = 0; col < N; col++) {
            for (int row = 0; row < N; row++) {
                // compute euclidean distance in dim dimensions
                squared_distance = 0;
                for (int d = 0; d < dim; d++) squared_distance += (trajectory1[d][col] - trajectory2[d][row]) * (trajectory1[d][col] - trajectory2[d][row]);

                // color recurrence points in black, non-recurrence points in white
                if(squared_distance <= squared_threshold) off_Image.setRGB(col, row, 0xFFFFFFFF);
                else off_Image.setRGB(col, row, 0xFF000000);

            }
        }
        return off_Image;
    }

    /**
     * @return the average vertical line length as an estimate of the minimal diagonal line length.
     */
    public int estimate_l_min(){
        long sum_v_line_lengths = 0;
        long num_v_lines = 0;
        for (int v = 0; v < v_hist.length; v++) {
            sum_v_line_lengths += v_hist[v] * v;
            num_v_lines += v_hist[v];
        }
        // should never be one.
        return (int) Math.max(2, Math.round(1.0 * sum_v_line_lengths / num_v_lines));
    }

    public void computeRQA(int minLengthDiagonal, int minLengthVertical, int minLengthOrthogonal){
        DRQA.computeRQA(minLengthDiagonal, minLengthVertical, minLengthOrthogonal,
                diagonal_rqa, vertical_rqa, white_vertical_rqa, orthogonal_rqa, diagonal_rqa_definite, vertical_rqa_definite, white_vertical_rqa_definite, orthogonal_rqa_definite,
                l_hist, l_hist_indefinite, v_hist, v_hist_indefinite, w_hist, w_hist_indefinite, r_hist, r_hist_indefinite);
    }
    public static void computeRQA(int minLengthDiagonal, int minLengthVertical, int minLengthOrthogonal,
                                  HashMap<HistogramStatistic, Double> diagonal_rqa, HashMap<HistogramStatistic, Double> vertical_rqa, HashMap<HistogramStatistic, Double> white_vertical_rqa, HashMap<HistogramStatistic, Double> orthogonal_rqa, HashMap<HistogramStatistic, Double> diagonal_rqa_definite, HashMap<HistogramStatistic, Double> vertical_rqa_definite, HashMap<HistogramStatistic, Double> white_vertical_rqa_definite, HashMap<HistogramStatistic, Double> orthogonal_rqa_definite,
                                  long[] l_hist, long[] l_hist_indefinite, long[] v_hist, long[] v_hist_indefinite, long[] w_hist, long[] w_hist_indefinite, long[] r_hist, long[] r_hist_indefinite){

        long[] histogram_definite, histogram_indefinite;
        int minLength;
        Map<HistogramStatistic, Double> resultMap;

        // compute for RQA and DRQA
        for(HistogramStatisticKind statisticsKind : HistogramStatisticKind.values()){

            // compute for all line types
            for(LineType lineType : LineType.values()){

                // depending on line type, the result map is changed (diagonal_rqa, vertical_rqa, etc.)
                switch (lineType){
                    case DIAGONAL:
                        histogram_definite = l_hist; histogram_indefinite = l_hist_indefinite;
                        minLength = minLengthDiagonal;
                        resultMap = statisticsKind == HistogramStatisticKind.RQA ? diagonal_rqa : diagonal_rqa_definite;
                        break;
                    case VERTICAL:
                        histogram_definite = v_hist; histogram_indefinite = v_hist_indefinite;
                        minLength = minLengthVertical;
                        resultMap = statisticsKind == HistogramStatisticKind.RQA ? vertical_rqa : vertical_rqa_definite;
                        break;
                    case WHITE_VERTICAL:
                        histogram_definite = w_hist; histogram_indefinite = w_hist_indefinite;
                        minLength = 1;
                        resultMap = statisticsKind == HistogramStatisticKind.RQA ? white_vertical_rqa : white_vertical_rqa_definite;
                        break;
                    default: // orthogonal lines
                        histogram_definite = r_hist; histogram_indefinite = r_hist_indefinite;
                        minLength = minLengthOrthogonal;
                        resultMap = statisticsKind == HistogramStatisticKind.RQA ? orthogonal_rqa : orthogonal_rqa_definite;
                }


                // Standard RQA doesn't differentiate between lines of definite and indefinite length.
                // To compute standard RQA measures, we simply sum the two histograms to level the partitioning between lines that are incident with the matrix border and lines that are not.
                long[] hist = statisticsKind == HistogramStatisticKind.RQA ? array_sum(histogram_definite, histogram_indefinite) : histogram_definite;
//                System.out.println(String.format("\nlineType: %s statisticsKind: %s", lineType, statisticsKind));
//                System.out.println(String.format("Arrays.toString(hist): %s", Arrays.toString(hist)));

                // basic aggregation -----------------------------------------------------------------------------------
                // number of lines larger than minimum line length
                long num_lines_filtered = 0; for (int i = minLength; i < hist.length; i++) num_lines_filtered += hist[i];
                // number of recurrence points on lines larger than minimum length
                long num_points_filtered = 0;
                for (int i = minLength; i < hist.length; i++) num_points_filtered += i * histogram_definite[i] + i * histogram_indefinite[i];
                // number of recurrence points on the lines in the histogram
                long num_points_unfiltered = 0;
                for (int i = 0; i < hist.length; i++) num_points_unfiltered += i * histogram_definite[i] + i * histogram_indefinite[i];

                // normalize sum of bins larger than minimum line length to one
                double[] filtered_normalized = new double[hist.length]; for (int i = minLength; i < hist.length; i++) filtered_normalized[i] = 1.*hist[i] / num_lines_filtered;
                // largest bin count
                long max_count = 0;
                for (long bin_count : hist) max_count = Math.max(bin_count, max_count);

                // -----------------------------------------------------------------------------------------------------
                // (D)RQA aggregation
                // -----------------------------------------------------------------------------------------------------

                // FILTER_RATIO:
                resultMap.put(HistogramStatistic.FILTER_RATIO, 1. * num_points_filtered / num_points_unfiltered);
//                System.out.println(String.format("num_points_filtered: %s", num_points_filtered));
//                System.out.println(String.format("num_points_unfiltered: %s", num_points_unfiltered));

                // AVERAGE: ∑ l * (P(l) / ∑ P(l))
                double average = 0; for (int i = minLength; i < filtered_normalized.length; i++) average += i * filtered_normalized[i];
                resultMap.put(HistogramStatistic.AVERAGE, Double.isNaN(average) ? 0 : average);

                // MEDIAN: empirical CDF(median) >= 0.5
                double median = 0;
                double mass_so_far = 0;
                for (int i = minLength; i < filtered_normalized.length; i++){
                    mass_so_far += filtered_normalized[i];
                    if(mass_so_far >= 0.5){
                        median = i;
                        break;
                    }
                }
                resultMap.put(HistogramStatistic.MEDIAN, Double.isNaN(median) ? 0 : median);

                // MAXIMUM: max{l | P(l) > 0}. By convention: zero if histogram is empty
                double max_line_length = 0;
                for (int i = hist.length - 1; i >= 0; i--) if(hist[i] > 0){ max_line_length = i; break; }
                resultMap.put(HistogramStatistic.MAX, Double.isNaN(max_line_length) ? 0 : max_line_length);

                // INVERSE_MAX or 0 if maximum line length is zero
                final double inverse_max = max_line_length > 0 ? 1. / max_line_length : 0;
                resultMap.put(HistogramStatistic.INVERSE_MAX, Double.isNaN(inverse_max) ? 0 : inverse_max);

                // ENTROPY: ∑ -p(l) ln p(l)
                double entropy = 0;
                for (int i = minLength; i < hist.length; i++) if(filtered_normalized[i] > 0) entropy += -1. * filtered_normalized[i] * Math.log(filtered_normalized[i]);
                resultMap.put(HistogramStatistic.ENTROPY, Double.isNaN(entropy) ? 0 : entropy);

                // DEFINITENESS: fraction points on lines of definite length
                long on_definite_lines = 0, on_indefinite_lines = 0;
                for (int i = 0; i < hist.length; i++) {
                    on_definite_lines += i * histogram_definite[i];
                    on_indefinite_lines += i * histogram_indefinite[i];
                }
                final double definiteness = 1. * on_definite_lines / (on_definite_lines + on_indefinite_lines);
                resultMap.put(HistogramStatistic.DEFINITENESS, Double.isNaN(definiteness) ? 0 : definiteness);

                // SORTEDNESS: intersection distance between sorted and unsorted histogram
                // sort only the part above minimum line length. Sort descending. It's faster to negate the elements and sort ascending (and negate at access) than using a DoubleComparator (which requires the array to be of boxed type Double[])
                double[] sorted = Arrays.copyOfRange(filtered_normalized, minLength, hist.length);
                for (int i = 0; i < sorted.length; i++) sorted[i] = -sorted[i];
                Arrays.sort(sorted);
                double sortedness = 0; for (int i = minLength; i < hist.length; i++) sortedness += Math.min(filtered_normalized[i], -sorted[i-minLength]);
                resultMap.put(HistogramStatistic.SORTEDNESS, Double.isNaN(sortedness) ? 0 : sortedness);

                // SPARSENESS: fraction of zero bins before maximum, zero if no bin can be zero (max line length = minimum line length)
                int zero_bins = 0;
                for (int i = minLength; i < max_line_length; i++)  if(hist[i] == 0) zero_bins++;
                final int possible_zero_bins = (int) max_line_length - minLength;
                final double sparseness_relative = 1. * zero_bins / (possible_zero_bins == 0 ? 1 : possible_zero_bins);
                resultMap.put(HistogramStatistic.SPARSENESS, Double.isNaN(sparseness_relative) ? 0 : sparseness_relative);

                // LOCAL_MAXIMA: number of local maxima in the histogram relative to the maximum number of local maxima (each second bin)
                int peaks = 0;
                int maxPossiblePeaks = (int) Math.ceil((float)(max_line_length-minLength+1)/2.0); // number of bins max-min+1 (b); each local maximum needs one subsequent element to be larger than (b/2); the last one is always followed by 0 (ceil)
                for (int i = minLength+1; i<hist.length-1; i++)  if(hist[i-1] < hist[i] && hist[i] > hist[i+1]) peaks++;
                if(hist[minLength] > hist[minLength+1]) peaks++;
                if(hist[hist.length-2] < hist[hist.length-1]) peaks++;  // since we assume that after the maximum line length, all bins are filled with zero, we can count the last element as local maximum if it is larger than its predecessor (its successor is zero and will be smaller than any non-zero element).
                final double local_maxima_relative = 1. * peaks / maxPossiblePeaks;
                resultMap.put(HistogramStatistic.LOCAL_MAXIMA, Double.isNaN(local_maxima_relative) ? 0 : local_maxima_relative);

                // DIFFERENCE ENTROPY: normalized entropy of the absolute difference between neighboring bins
                SparseDoubleMatrix1D stepCounts = new SparseDoubleMatrix1D(hist.length);
                int stepsTotal = (int) max_line_length - minLength;
                long max_step_size = 0;
                for (int i = minLength+1; i <= max_line_length; i++) {
                    int step = (int) Math.abs(hist[i] - hist[i-1]);
                    max_step_size = step > max_step_size ? step : max_step_size;
                    double previousCount = stepCounts.getQuick(step);
                    stepCounts.setQuick(step, previousCount+1./stepsTotal);
                }
                double diff_entr = 0;
                IntArrayList indexList = new IntArrayList(); DoubleArrayList valueList = new DoubleArrayList();
                stepCounts.getNonZeros(indexList, valueList);
                long largestDifference = 0;
                for (int i = 0; i < valueList.size(); i++){
                    largestDifference = Math.max(largestDifference, (long) valueList.get(i));
                    diff_entr += - valueList.get(i) * Math.log(valueList.get(i));
                }
                if(largestDifference > 1) diff_entr /= Math.log(largestDifference);
                resultMap.put(HistogramStatistic.DIFFERENCE_ENTROPY, Double.isNaN(diff_entr) ? 0 : diff_entr);

            }
        }

    }

    /**
     * Computes the maximum euclidean distance between any point on the first and any point on the second trajectory.
     * Quadratic runtime.
     * @param trajectory1 First trajectory. The array dimensions correspond to [phase space dimension][time step].
     * @param trajectory2 Second trajectory.
     * @return The maximum distance between any two points.
     */
    public static double maxPhaseSpaceDiameter(double[][] trajectory1, double[][] trajectory2){
        final int N = trajectory1[0].length;
        double max_squared_distance = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {
                double squared_distance = 0;
                for (int d = 0; d < trajectory1.length; d++) squared_distance += (trajectory1[d][i] - trajectory2[d][j]) * (trajectory1[d][i] - trajectory2[d][j]);
                max_squared_distance = Math.max(max_squared_distance, squared_distance);
            }
        }
        return Math.sqrt(max_squared_distance);
    }

    /**
     * Approximates the phase space diamter by computing the diagonal of the bounding (axis parallel) hypercube.
     * @param ts The input time series [phase space dimension][time step]
     * @return An upper bound of the maximum phase space diameter.
     */
    public static double maxPhaseSpaceDiameterApproximate(double[][] ts){

        int numDims = ts.length;
        double[] mins = new double[numDims], maxs = new double[numDims];
        Arrays.fill(mins, Double.POSITIVE_INFINITY);
        Arrays.fill(maxs, Double.NEGATIVE_INFINITY);

        // find the minimum and maximum in each dimension
//        double[] mean = new double[numDims];
        for (int dim = 0; dim < numDims; dim++) {
            for (int t = 0; t < ts[dim].length; t++) {
                if(ts[dim][t] < mins[dim]) mins[dim] = ts[dim][t];
                if(ts[dim][t] > maxs[dim]) maxs[dim] = ts[dim][t];
//                mean[dim] += ts[dim][t];
            }
//            mean[dim] /= ts[dim].length;
        }


        // The length of all diagonals of the hypercube is the same (because one needs to make the step in each dimension to reach the opposite point).
        // The square root of the sum of squared side lengths.
        double sum = 0;
        for (int dim = 0; dim < numDims; dim++) {
            sum += (maxs[dim] - mins[dim]) * (maxs[dim] - mins[dim]);
        }
        return Math.sqrt(sum);

    }

    public double getDET_RR_RATIO(){
        double det_rr_ratio = diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO) / recurrence_rate;
        return Double.isNaN(det_rr_ratio) ? -1 : det_rr_ratio;
    }

    public double getLAM_DET_RATIO(){
        double lam_det_ratio = vertical_rqa_definite.get(HistogramStatistic.FILTER_RATIO) / diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO);
        return Double.isNaN(lam_det_ratio) ? -1 : lam_det_ratio;
    }

    /**
     * Used as a fast replacement for histogram description in bootstrapping.
     */
    public static void computeRQAForHist(int minLength, final long[] histogram_definite, double[] result){

        // basic aggregation -----------------------------------------------------------------------------------
        // number of lines larger than minimum line length
        long num_lines_filtered = 0; for (int i = minLength; i < histogram_definite.length; i++) num_lines_filtered += histogram_definite[i];
        // number of recurrence points on lines larger than minimum length
        long num_points_filtered = 0;
        for (int i = minLength; i < histogram_definite.length; i++) num_points_filtered += i * histogram_definite[i];
        // number of recurrence points on the lines in the histogram
        long num_points_unfiltered = 0;
        for (int i = 0; i < histogram_definite.length; i++) num_points_unfiltered += i * histogram_definite[i];
        // normalize sum of bins larger than minimum line length to one
        double[] filtered_normalized = new double[histogram_definite.length]; for (int i = minLength; i < histogram_definite.length; i++) filtered_normalized[i] = 1.*histogram_definite[i] / num_lines_filtered;

        // -----------------------------------------------------------------------------------------------------
        // RQA aggregation
        // -----------------------------------------------------------------------------------------------------

        // FILTER_RATIO:
        double filter_ratio = 1. * num_points_filtered / num_points_unfiltered;

        // AVERAGE: ∑ l * (P(l) / ∑ P(l))
        double average = 0; for (int i = minLength; i < filtered_normalized.length; i++) average += i * filtered_normalized[i];

        // MEDIAN: empirical CDF(median) >= 0.5
        double median = 0;
        double mass_so_far = 0;
        for (int i = minLength; i < filtered_normalized.length; i++){
            mass_so_far += filtered_normalized[i];
            if(mass_so_far >= 0.5){
                median = i;
                break;
            }
        }

        // ENTROPY: ∑ -p(l) ln p(l)
        double entropy = 0;
        for (int i = minLength; i < histogram_definite.length; i++) if(filtered_normalized[i] > 0) entropy += -1. * filtered_normalized[i] * Math.log(filtered_normalized[i]);

        result[0] = filter_ratio;
        result[1] = average;
        result[2] = median;
        result[3] = entropy;

    }

    /**
     * Estimates the mean and confidence intervals of filter ratio, average, median and entropy for a given line length histogram.
     * @param p_value The confidence interval p-value.
     * @param resamples The number of line length histograms to generate by bootstrapping the given histogram. 1000 is a common value, although for complex distributions 10000 returned more stable results (takes long for large histograms).
     * @param min_line_length The minimum line length for the filter ratio, average, median and entropy.
     * @param line_length_distribution The histogram to use as base for bootstrapping.
     * @param mean A double array of size four. Used to store the mean of the sampling distributions of filter ratio, average, median, and entropy.
     * @param ci A double array of size four. Used to store the confidence interval width. This is the difference between the (1-p/2)-percentile and the (p/2)-percentile of the sampling distributions of filter ratio, average, median, and entropy.
     */
    public static void bootstrapSamplingDistribution(double p_value, int resamples, int min_line_length, long[] line_length_distribution, double[] mean, double[] ci){

        if(resamples == 0) return;

        // to summarize the sampling distribution of filter ratio, average, median and entropy
        DescriptiveStatistics[] sampling_distributions = new DescriptiveStatistics[4];
        for (int i = 0; i < sampling_distributions.length; i++)  sampling_distributions[i] = new DescriptiveStatistics(resamples+50);

        // the total number of lines recorded in the histogram
        long num_lines = 0; for (long bin_value : line_length_distribution) num_lines += bin_value;
        // the normalized histogram (all bins sum to one)
        double[] normalized = new double[line_length_distribution.length]; for (int i = 0; i < line_length_distribution.length; i++) normalized[i] = 1.*line_length_distribution[i] / num_lines;

        // the empirical probability distribution from which "lines" are sampled
        Empirical empirical = new Empirical(normalized, Empirical.NO_INTERPOLATION, new DRand());

        // contains the filter ratio, average, median and entropy for a single bootstrap sample
        double[] measures_temp = new double[4];
        long[] bootstrapped = new long[normalized.length];
        for (int i = 0; i < resamples; i++) {
//            if(i%10000==0) System.out.println(i);
            // generate a bootstrap sample
            Arrays.fill(bootstrapped, 0);
            for (int line_idx = 0; line_idx < num_lines; line_idx++) {
                int random_line_length = (int)(empirical.nextDouble()*normalized.length);
                bootstrapped[random_line_length]++;
            }
            // compute rqa and record results
            computeRQAForHist(min_line_length, bootstrapped, measures_temp);
            for (int j = 0; j < measures_temp.length; j++) sampling_distributions[j].addValue(measures_temp[j]);
        }

        // report mean and confidence interval width for each sampling distribution
        for (int i = 0; i < measures_temp.length; i++) {
            mean[i] = sampling_distributions[i].getMean();
            ci[i] = sampling_distributions[i].getPercentile(1.0 - p_value/2.0) - sampling_distributions[i].getPercentile(p_value/2.0);
        }

    }

    public static String[] allRQAMeasureNames(){

        String[] allNames = new String[LineType.values().length * HistogramStatistic.values().length];

        Map<LineType, Map<HistogramStatistic, Double>> maps = new HashMap<LineType, Map<HistogramStatistic, Double>>();

        int index = 0;
        for(LineType lineType : LineType.values())
            for(HistogramStatistic statistic : HistogramStatistic.values()){
                allNames[index] = lineType.toString() + "_" + statistic.toString();
                index++;
            }

        return allNames;
    }

    public double[] standardRQAMeasures(){
        return new double[]{
            recurrence_rate,
            diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO),                     // DET
            diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO)/recurrence_rate,     // DET/RR
            diagonal_rqa_definite.get(HistogramStatistic.AVERAGE),
            diagonal_rqa_definite.get(HistogramStatistic.MAX),
            diagonal_rqa_definite.get(HistogramStatistic.INVERSE_MAX),
            diagonal_rqa_definite.get(HistogramStatistic.ENTROPY),

            vertical_rqa_definite.get(HistogramStatistic.FILTER_RATIO),
            vertical_rqa_definite.get(HistogramStatistic.FILTER_RATIO)/diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO),  // LAM/DET
            vertical_rqa_definite.get(HistogramStatistic.AVERAGE),
            vertical_rqa_definite.get(HistogramStatistic.MAX),
            vertical_rqa_definite.get(HistogramStatistic.ENTROPY),

            white_vertical_rqa_definite.get(HistogramStatistic.AVERAGE),
            white_vertical_rqa_definite.get(HistogramStatistic.MAX),
            white_vertical_rqa_definite.get(HistogramStatistic.INVERSE_MAX),
            white_vertical_rqa_definite.get(HistogramStatistic.ENTROPY)
        };
    }
    /**
     * @param selection the selection to include in the result array
     * @return all rqa measures in a double array. The order is (l, v, w, r) x (FILTER_RATIO, AVERAGE,  MAX, INVERSE_MAX, ENTROPY, DEFINITENESS, SORTEDNESS, LOCAL_MAXIMA, SPARSENESS, DIFFERENCE_ENTROPY)
     */
    public double[] allRQAMeasures(HistogramStatistic[] selection){

        if(selection == null) selection = HistogramStatistic.values();

        double[] allMeasures = new double[LineType.values().length * selection.length];

        Map<LineType, Map<HistogramStatistic, Double>> maps = new HashMap<LineType, Map<HistogramStatistic, Double>>();
        maps.put(LineType.DIAGONAL, diagonal_rqa_definite);
        maps.put(LineType.VERTICAL, vertical_rqa_definite);
        maps.put(LineType.WHITE_VERTICAL, white_vertical_rqa_definite);
        maps.put(LineType.ORTHOGONAL, orthogonal_rqa_definite);

        int index = 0;
        for(LineType lineType : LineType.values()){
            Map<HistogramStatistic, Double> map = maps.get(lineType);
            for(HistogramStatistic statistic : selection){
                allMeasures[index] = map.get(statistic);
                index++;
            }
        }

        return allMeasures;

    }

    private static long[] array_sum(long[] array1, long[] array2) {
        final long[] result = new long[array1.length];
        for (int i = 0; i < array1.length; i++) result[i] = array1[i] + array2[i];
        return result;
    }

    /**
     * Divides each bin (an absolute count) by the total sum of counts (results in a frequency, all bins sum to one).
     * @param hist Counts
     * @return Frequencies
     */
    private double[] array_normalize(long[] hist){
        long sum = 0;
        for (long countInBin : hist) sum += countInBin;
        double[] result = new double[hist.length];
        for (int i = 0; i < hist.length; i++)  result[i] = 1. * hist[i] / sum;
        return result;
    }

    /**
     * Checks whether the sum of recurrence points (recurrence mass) matches the individually counted number of recurrence points.
     * Enable assertions of the JVM (-ea) to be notified in case of some deviation.
     * @param N the length of the time series / number of rows and columns in the RP.
     */
    private void checkHistogramSums(int N) {

        // the diagonal lines are counted only once, so to compare checksums, we count each line twice (except the main diagonal)
        long recurrence_points_comparison = 0;
        for (int l = 0; l < l_hist.length; l++) recurrence_points_comparison += 2 * l * l_hist[l] + 2 * l * l_hist_indefinite[l];
        assert recurrence_points == recurrence_points_comparison - N: String.format("Diagonal lines don't sum up. Recurrence points: %s sum of diagonal lines: %s", recurrence_points, recurrence_points_comparison);

        // the vertical line should sum up to the number of recurrence points
        recurrence_points_comparison = 0;
        for (int v = 0; v < v_hist.length; v++) recurrence_points_comparison += v * v_hist[v] + v * v_hist_indefinite[v];
        assert recurrence_points == recurrence_points_comparison: String.format("Vertical lines don't sum up. Recurrence points: %s sum of vertical lines: %s", recurrence_points, recurrence_points_comparison);

        // the white lines should sum up to the number of non-recurrence points (N*N - recurrence points)
        recurrence_points_comparison = 0;
        for (int w = 0; w < w_hist.length; w++) recurrence_points_comparison += w * w_hist[w] + w * w_hist_indefinite[w];
        assert N*N - recurrence_points == recurrence_points_comparison: String.format("White vertical lines don't sum up. Recurrence points: %s sum of white vertical lines: %s", recurrence_points, recurrence_points_comparison);

        // as for diagonal lines, we count each orthogonal line twice. Since the main diagonal is excluded, we add these N points (they miss either as trivial orthogonal lines of length one, or as missing points on main-diagonal-adjacent orthogonal lines).
        recurrence_points_comparison = 0;
        for (int r = 0; r < r_hist.length; r++) recurrence_points_comparison += 2 * r * r_hist[r] + 2 * r * r_hist_indefinite[r];
        assert recurrence_points == recurrence_points_comparison + N: String.format("Orthogonal lines don't sum up. Recurrence points: %s sum of orthogonal lines: %s", recurrence_points, recurrence_points_comparison);
    }

    public StringBuilder conditionalRecurrenceTsv(SparseDoubleMatrix2D conditional_recurrence, int limit){
        // write out conditional recurrence distribution
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < limit; i++) {
            for (int j = 0; j < limit; j++)
                builder.append(String.format("%s\t", (long) conditional_recurrence.getQuick(i, j)));
            builder.append("\n");
        }
        return builder;
    }

    public void writeCRTImages(String outpath){

        SparseDoubleMatrix2D[] conditional_recurrences = new SparseDoubleMatrix2D[]{conditional_ww}; // , conditional_vv, conditional_wv};
        int[] limits = new int[]{conditional_ww_limit}; //, conditional_vv_limit, conditional_wv_limit};
        String[] filename_suffices = new String[]{"ww", "vv", "wv"};

        for (int cr_idx = 0; cr_idx < conditional_recurrences.length; cr_idx++) {

            if(limits[cr_idx] <= 0) continue;

            SparseDoubleMatrix2D conditional_recurrence = conditional_recurrences[cr_idx];
            BufferedImage image_out = getCRTImage(limits[cr_idx], conditional_recurrence);

            String[] tokens = outpath.split("\\.(?=[^\\.]+$)");
            final String filename = tokens[0] + filename_suffices[cr_idx] + ".png";
            File outputfile = new File(filename);
            try {
                ImageIO.write(image_out, "png", outputfile);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public BufferedImage getCRTImage(int limit, SparseDoubleMatrix2D conditional_recurrence) {
        BufferedImage image_out = new BufferedImage(limit, limit, BufferedImage.TYPE_INT_ARGB);

        IntArrayList indexXList = new IntArrayList();
        IntArrayList indexYList = new IntArrayList();
        DoubleArrayList valueList = new DoubleArrayList();

        // find maximum to scale to the color range of the image
        conditional_recurrence.getNonZeros(indexXList, indexYList, valueList);
        double max = Double.NEGATIVE_INFINITY;
        for (double d : valueList.elements()) max = Math.max(d, max);

        if (CRT_LOG_SCALE) max = Math.log10(max);

        for (int i = 0; i < limit; i++) {
            for (int j = 0; j < limit; j++) {
                double value = conditional_recurrence.getQuick(i, j);
                if (CRT_LOG_SCALE) value = Math.log10(value);
                float color = (float) (value / max);
                image_out.setRGB(i, j, Color.HSBtoRGB(0f, 0f, color));
            }
        }
        return image_out;
    }

//     * Returns the CRT image but subtracts the image that results from the independence assumption.
//     * Since the observed number might be smaller than the expected number, color is used in this image.
//     * Black corresponds to zero, blue corresponds to negative numbers and red to positive numbers.
//     * @param limit
//     * @param conditional_recurrence
//     * @return
//     */
//    public BufferedImage getCRTImageDifferenceToIndependent(int limit, SparseDoubleMatrix2D conditional_recurrence) {
//        BufferedImage image_out = new BufferedImage(limit, limit, BufferedImage.TYPE_INT_ARGB);
//
//        double[][] crt = conditional_recurrence.toArray();
//        // compute the difference histogram and find maximum to scale to the color range
//        double max = Double.NEGATIVE_INFINITY;
//        int rows = crt.length;
//        int cols = crt[0].length;
//        int numLines = conditional_recurrence.cardinality();
//        for (int i = 0; i < rows; i++) {
//            for (int j = 0; j < cols; j++) {
//                crt[i][j] -=
//                max = Math.max(Math.abs(d), max);
//            }
//        }
//
//        if (CRT_LOG_SCALE) max = Math.log10(max);
//
//        for (int i = 0; i < limit; i++) {
//            for (int j = 0; j < limit; j++) {
//                double value = conditional_recurrence.getQuick(i, j);
//                if (CRT_LOG_SCALE) value = Math.log10(value);
//                float brightness = (float) (value / max);
//
//                image_out.setRGB(i, j, Color.HSBtoRGB(0f, 0f, brightness));
//            }
//        }
//        return image_out;
//    }
//    /**
//     * Returns the CRT image but subtracts the image that results from the independence assumption.
//     * Since the observed number might be smaller than the expected number, color is used in this image.
//     * Black corresponds to zero, blue corresponds to negative numbers and red to positive numbers.
//     * @param limit
//     * @param conditional_recurrence
//     * @return
//     */
//    public BufferedImage getCRTImageDifferenceToIndependent(int limit, SparseDoubleMatrix2D conditional_recurrence) {
//        BufferedImage image_out = new BufferedImage(limit, limit, BufferedImage.TYPE_INT_ARGB);
//
//        double[][] crt = conditional_recurrence.toArray();
//        // compute the difference histogram and find maximum to scale to the color range
//        double max = Double.NEGATIVE_INFINITY;
//        int rows = crt.length;
//        int cols = crt[0].length;
//        int numLines = conditional_recurrence.cardinality();
//        for (int i = 0; i < rows; i++) {
//            for (int j = 0; j < cols; j++) {
//                crt[i][j] -=
//                max = Math.max(Math.abs(d), max);
//            }
//        }
//
//        if (CRT_LOG_SCALE) max = Math.log10(max);
//
//        for (int i = 0; i < limit; i++) {
//            for (int j = 0; j < limit; j++) {
//                double value = conditional_recurrence.getQuick(i, j);
//                if (CRT_LOG_SCALE) value = Math.log10(value);
//                float brightness = (float) (value / max);
//
//                image_out.setRGB(i, j, Color.HSBtoRGB(0f, 0f, brightness));
//            }
//        }
//        return image_out;
//    }


    public void print(){

        Locale.setDefault(Locale.US);
        System.out.println(String.format("recurrence_rate: %s", recurrence_rate));
        Map<HistogramStatistic, Double>[] maps = new Map[]{diagonal_rqa_definite, vertical_rqa_definite, white_vertical_rqa_definite, orthogonal_rqa_definite};
        String[] names = new String[]{"diagonal_rqa_definite", "vertical_rqa_definite", "white_vertical_rqa_definite", "orthogonal_rqa_definite"};

        for (int i = 0; i < maps.length; i++) {
            System.out.println(names[i]);

            System.out.println(String.format("FILTER_RATIO......... %.3f", maps[i].get(HistogramStatistic.FILTER_RATIO)));
            System.out.println(String.format("AVERAGE.............. %.3f", maps[i].get(HistogramStatistic.AVERAGE)));
            System.out.println(String.format("MEDIAN............... %s", Math.round(maps[i].get(HistogramStatistic.MEDIAN))));
            System.out.println(String.format("MAX.................. %s", Math.round(maps[i].get(HistogramStatistic.MAX))));
            System.out.println(String.format("ENTROPY.............. %.3f", maps[i].get(HistogramStatistic.ENTROPY)));
            System.out.println(String.format("SORTEDNESS........... %.3f", maps[i].get(HistogramStatistic.SORTEDNESS)));
            System.out.println(String.format("DIFFERENCE.ENTROPY... %.3f", maps[i].get(HistogramStatistic.DIFFERENCE_ENTROPY)));
            System.out.println(String.format("SPARSENESS........... %.3f", maps[i].get(HistogramStatistic.SPARSENESS)));
            System.out.println(String.format("LOCAL_MAXIMA......... %.3f", maps[i].get(HistogramStatistic.LOCAL_MAXIMA)));
            System.out.println(String.format("DEFINITENESS......... %.3f\n", maps[i].get(HistogramStatistic.DEFINITENESS)));
        }

    }

    public void printComparison(DRQA other) {

        Locale.setDefault(Locale.US);
        System.out.println(String.format("recurrence_rate: %s", recurrence_rate));
        Map<HistogramStatistic, Double>[] maps = new Map[]{diagonal_rqa_definite, vertical_rqa_definite, white_vertical_rqa_definite, orthogonal_rqa_definite};
        Map<HistogramStatistic, Double>[] other_maps = new Map[]{other.diagonal_rqa_definite, other.vertical_rqa_definite, other.white_vertical_rqa_definite, other.orthogonal_rqa_definite,};
        String[] names = new String[]{"diagonal_rqa_definite", "vertical_rqa_definite", "white_vertical_rqa_definite", "orthogonal_rqa_definite"};

        for (int i = 0; i < maps.length; i++) {
            System.out.println(names[i]);

            System.out.println(String.format("FILTER_RATIO......... %.3f\t%.3f", maps[i].get(HistogramStatistic.FILTER_RATIO), other_maps[i].get(HistogramStatistic.FILTER_RATIO)));
            System.out.println(String.format("AVERAGE.............. %.3f\t%.3f", maps[i].get(HistogramStatistic.AVERAGE), other_maps[i].get(HistogramStatistic.AVERAGE)));
            System.out.println(String.format("MEDIAN............... %s\t%s", Math.round(maps[i].get(HistogramStatistic.MEDIAN)), Math.round(other_maps[i].get(HistogramStatistic.MEDIAN))));
            System.out.println(String.format("MAX.................. %s\t%s", Math.round(maps[i].get(HistogramStatistic.MAX)), Math.round(other_maps[i].get(HistogramStatistic.MAX))));
            System.out.println(String.format("ENTROPY.............. %.3f\t%.3f", maps[i].get(HistogramStatistic.ENTROPY), other_maps[i].get(HistogramStatistic.ENTROPY)));
            System.out.println(String.format("SORTEDNESS........... %.3f\t%.3f", maps[i].get(HistogramStatistic.SORTEDNESS), other_maps[i].get(HistogramStatistic.SORTEDNESS)));
            System.out.println(String.format("DIFFERENCE.ENTROPY... %.3f\t%.3f", maps[i].get(HistogramStatistic.DIFFERENCE_ENTROPY), other_maps[i].get(HistogramStatistic.DIFFERENCE_ENTROPY)));
            System.out.println(String.format("SPARSENESS........... %.3f\t%.3f", maps[i].get(HistogramStatistic.SPARSENESS), other_maps[i].get(HistogramStatistic.SPARSENESS)));
            System.out.println(String.format("LOCAL_MAXIMA......... %.3f\t%.3f", maps[i].get(HistogramStatistic.LOCAL_MAXIMA), other_maps[i].get(HistogramStatistic.LOCAL_MAXIMA)));
            System.out.println(String.format("DEFINITENESS......... %.3f\t%.3f\n", maps[i].get(HistogramStatistic.DEFINITENESS), other_maps[i].get(HistogramStatistic.DEFINITENESS)));
        }
    }

    public static String histogramToString(long[] v_hist) {
        int lastNonZeroBinIdx = v_hist.length-1;
        for (int i = v_hist.length-1; i >= 0; i--) {
            if(v_hist[i] > 0) break;
            lastNonZeroBinIdx = i;
        }

        return Arrays.toString(Arrays.copyOf(v_hist,lastNonZeroBinIdx)).replaceAll(" |\\[|\\]","");
    }

    /** Is used in text output to separate columns. */
    public static String separator = "\t";
    public static String separatedHeader(){
        final String measureNames = Arrays.toString(DRQA.allRQAMeasureNames()).replaceAll(" |\\[|\\]", "").replaceAll(",", separator);
        return String.format("path%1$sclass_name%1$sembedding_method_id%1$spst_dimension%1$sembedding_delay%1$strajectory_length%1$s" +
                "recurrence_threshold%1$snorm%1$sL_min_method_id%1$sL_min%1$snoise_ratio_id%1$s" +
                "recurrence_rate%1$sDET_RR_RATIO%1$sLAM_DET_RATIO%1$s", separator) + measureNames + String.format("%1$sL_HIST%1$sV_HIST%1$sW_HIST%1$sR_HIST%1$scrt_mean_row%1$scrt_mean_col%1$scrt_correlation%1$scrt_max_row%1$scrt_max_col%1$scrt_local_maxima%1$scrt_entropy", separator);

    }

    /**
     * @param filename The identifier of the underlying data.
     * @param embeddingMethodId The embedding method, e.g. original data, no embedding or time delay embedded.
     * @param dimension The dimensionality of the phase space trajectory.
     * @param delay The delay in embedding. Should be -1 if no embedding was applied to create the phase space trajectory.
     * @param l_min minimum diagonal line length.
     * @param noiseRatioId The noise ratio that was applied to the trajectory.
     * @param trajectory_length the length of the trajectory
     * @return Creates a single, separated string for output in a separated text file.
     */
    public String toSeparatedString(String filename, String class_name, int embeddingMethodId, int dimension, int delay, int l_min_method_id, int l_min, int noiseRatioId, int trajectory_length) {
        double[] allMeasures = allRQAMeasures(HistogramStatistic.values());
        String measures_tsv = Arrays.toString(allMeasures).replaceAll("\\[|\\]","").replaceAll(",", separator);
        String histograms_tsv =
                histogramToString(l_hist) + separator +
                        histogramToString(v_hist) + separator +
                        histogramToString(w_hist) + separator +
                        histogramToString(r_hist);
        return String.format("%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s%18$s%s",
                filename, class_name, embeddingMethodId, dimension, delay, trajectory_length, recurrence_threshold, "euclidean", l_min_method_id, l_min, noiseRatioId,
                recurrence_rate, getDET_RR_RATIO(), getLAM_DET_RATIO(), measures_tsv, histograms_tsv, crtStatistics(), separator);
    }

    public enum EmbeddingMethod{
        ORIGINAL_TRAJECTORY, TIME_DELAY_EMBEDDING, ONE_DIMENSION, TIME_DELAY_EMBEDDING_DIM2, TIME_DELAY_EMBEDDING_SMALLEST_MIN, TIME_DELAY_EMBEDDING_FIRST_MIN_UNBOUNDED
    }

    public enum LMinMethod{
        L_MIN_FIX, L_MIN_EST
    }

    public void writeResultDB(String path, String class_name, int embeddingMethodId, int dimension, int delay, int trajectory_length, int l_min_method_id, int l_min, int noiseRatioId, Connection connection) {
        String insert_sql = "INSERT INTO `db_wittcarl`.`rqa_measures`" +
                // 1           2                  3                   4               5                     6                    7                 8              9           10          11           12          13                14                 15                  16                  17                    18                 19                  20                  21                      22                        23                      24                       25                          26                          27                28             29                      30                31                   32                    33                        34                          35                      36                               37                         38                        39                       40                      41                      42                      43                                44                               45                              46                          47                                 48                       49                     50                  51                    52                   53                        54                          55                      56                      57                          58                      59            60      61        62      63
                "(`path`, `class_name`, `embedding_method_id`, `pst_dimension`, `embedding_delay`, `trajectory_length`, `recurrence_threshold`, `norm`, `l_min_method_id`, `l_min`, `noise_ratio_id`, `RR`, `DET_RR_RATIO`, `LAM_DET_RATIO`, `DIAGONAL_FILTER_RATIO`, `DIAGONAL_AVERAGE`, `DIAGONAL_MEDIAN`, `DIAGONAL_MAX`, `DIAGONAL_INVERSE_MAX`, `DIAGONAL_ENTROPY`, `DIAGONAL_DEFINITENESS`, `DIAGONAL_SORTEDNESS`, `DIAGONAL_LOCAL_MAXIMA`, `DIAGONAL_SPARSENESS`, `DIAGONAL_DIFFERENCE_ENTROPY`, `VERTICAL_FILTER_RATIO`, `VERTICAL_AVERAGE`, `VERTICAL_MEDIAN`, `VERTICAL_MAX`, `VERTICAL_INVERSE_MAX`, `VERTICAL_ENTROPY`, `VERTICAL_DEFINITENESS`, `VERTICAL_SORTEDNESS`, `VERTICAL_LOCAL_MAXIMA`, `VERTICAL_SPARSENESS`, `VERTICAL_DIFFERENCE_ENTROPY`, `WHITE_VERTICAL_FILTER_RATIO`, `WHITE_VERTICAL_AVERAGE`, `WHITE_VERTICAL_MEDIAN`, `WHITE_VERTICAL_MAX`, `WHITE_VERTICAL_INVERSE_MAX`, `WHITE_VERTICAL_ENTROPY`, `WHITE_VERTICAL_DEFINITENESS`, `WHITE_VERTICAL_SORTEDNESS`, `WHITE_VERTICAL_LOCAL_MAXIMA`, `WHITE_VERTICAL_SPARSENESS`, `WHITE_VERTICAL_DIFFERENCE_ENTROPY`, `ORTHOGONAL_FILTER_RATIO`, `ORTHOGONAL_AVERAGE`, `ORTHOGONAL_MEDIAN`, `ORTHOGONAL_MAX`, `ORTHOGONAL_INVERSE_MAX`, `ORTHOGONAL_ENTROPY`, `ORTHOGONAL_DEFINITENESS`, `ORTHOGONAL_SORTEDNESS`, `ORTHOGONAL_LOCAL_MAXIMA`, `ORTHOGONAL_SPARSENESS`, `ORTHOGONAL_DIFFERENCE_ENTROPY`, `L_HIST`, `V_HIST`, `W_HIST`, `R_HIST`, `signature`)" +
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, COMPRESS(?), COMPRESS(?), COMPRESS(?), ?)";
        PreparedStatement preps = null;

        try {
            preps = connection.prepareStatement(insert_sql);
            preps.setString(1, path);
            preps.setString(2, class_name);
            preps.setInt(3, embeddingMethodId);
            preps.setInt(4, dimension);
            preps.setInt(5, delay);
            preps.setInt(6, trajectory_length);
            preps.setFloat(7, ((float) recurrence_threshold));
            preps.setString(8, "euclidean");
            preps.setInt(9, l_min_method_id);
            preps.setInt(10, l_min);
            preps.setInt(11, noiseRatioId);
            preps.setFloat(12, ((float) (Double.isNaN(recurrence_rate) ? -1 : recurrence_rate)));
            final double det_rr_ratio = diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO) / recurrence_rate;
            preps.setFloat(13, (float) (Double.isNaN(det_rr_ratio) ? -1 : det_rr_ratio));
            final double lam_det_ratio = vertical_rqa_definite.get(HistogramStatistic.FILTER_RATIO) / diagonal_rqa_definite.get(HistogramStatistic.FILTER_RATIO);
            preps.setFloat(14, (float) (Double.isNaN(lam_det_ratio) ? -1 : lam_det_ratio));
            int nextIdx = 15;
            double[] rqa_measures = allRQAMeasures(null);
            for (int offset = 0; offset < rqa_measures.length; offset++)
                preps.setFloat(nextIdx+offset, (float) (Double.isNaN(rqa_measures[offset]) ? -1 : rqa_measures[offset]));
            preps.setString(59, histogramToString(l_hist));
            preps.setString(60, histogramToString(v_hist));
            preps.setString(61, histogramToString(w_hist));
            preps.setString(62, histogramToString(r_hist));
            preps.setObject(63, null);
            preps.execute();
        } catch (SQLException e) {
            e.printStackTrace();
        }


    }

    protected static short[][] compress_hist(SparseDoubleMatrix2D matrix) {
        double maxVal = Double.NEGATIVE_INFINITY;
        int maxRow = 0, maxCol = 0;
        IntArrayList rowList = new IntArrayList(), colList = new IntArrayList();
        DoubleArrayList valueList = new DoubleArrayList();
        matrix.getNonZeros(rowList, colList, valueList);
        for (int i = 0; i < valueList.size(); i++) maxVal = Math.max(maxVal, valueList.get(i));

        double maxValLog = Math.log10(maxVal);
        for (int i = 0; i < valueList.size(); i++) {
            if( Math.log10(valueList.get(i) + 1) > maxValLog / 20.0){
                maxRow = Math.max(maxRow, rowList.get(i));
                maxCol = Math.max(maxCol, colList.get(i));
            }
        }

        short[][] result = new short[maxRow+1][maxCol+1];
        for (int row = 0; row <= maxRow; row++) {
            for (int col = 0; col <= maxCol; col++) {
                result[row][col] = (short) Math.log10(matrix.getQuick(row, col)+1);
            }
        }

        return result;
    }

}

//c("LENGTH","L_FILTER_RATIO", "L_AVERAGE",  "L_MAX", "L_INVERSE_MAX", "L_ENTROPY", "L_DEFINITENESS", "L_SORTEDNESS", "L_LOCAL_MAXIMA", "L_SPARSENESS", "L_DIFFERENCE_ENTROPY",
//"V_FILTER_RATIO", "V_AVERAGE",  "V_MAX", "V_INVERSE_MAX", "V_ENTROPY", "V_DEFINITENESS", "V_SORTEDNESS", "V_LOCAL_MAXIMA", "V_SPARSENESS", "V_DIFFERENCE_ENTROPY",
//        "W_FILTER_RATIO", "W_AVERAGE",  "W_MAX", "W_INVERSE_MAX", "W_ENTROPY", "W_DEFINITENESS", "W_SORTEDNESS", "W_LOCAL_MAXIMA", "W_SPARSENESS", "W_DIFFERENCE_ENTROPY",
//        "R_FILTER_RATIO", "R_AVERAGE",  "R_MAX", "R_INVERSE_MAX", "R_ENTROPY", "R_DEFINITENESS", "R_SORTEDNESS", "R_LOCAL_MAXIMA", "R_SPARSENESS", "R_DIFFERENCE_ENTROPY")