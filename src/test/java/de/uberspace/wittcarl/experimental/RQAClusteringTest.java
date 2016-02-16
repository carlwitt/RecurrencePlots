package de.uberspace.wittcarl.experimental;

import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import junit.framework.TestCase;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.clustering.evaluation.ClusterEvaluator;
import org.apache.commons.math3.ml.clustering.evaluation.SumOfClusterVariances;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.random.ISAACRandom;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Author: Carl Witt, de.uberspace.wittcarl@deneb.uberspace.de
 * Date: 20.08.15.
 */
public class RQAClusteringTest extends TestCase{

    public void testClustering(){

        final int numberOfTimeSeriesPerClass = 100;
        final int timeSeriesLength = 1500;
        final int sample_size = (9 * numberOfTimeSeriesPerClass) / 2;
        final int numberOfResamples = 100;
        final int maxIterations = 1000;
        final int numberOfInitializations = 100;

        // ------------------------------------
        // compute DRQA measures
        // ------------------------------------

        TimeSeriesGenerator timeSeriesGenerator = new TimeSeriesGenerator(numberOfTimeSeriesPerClass, timeSeriesLength, null);
        ArrayList<DoublePoint> rqa_measures = new ArrayList<DoublePoint>(9*numberOfTimeSeriesPerClass);

        System.out.println("Computing DRQA measures");
        double[][][][] trajectories = timeSeriesGenerator.getAllTrajectories();
        for (int system_idx = 0; system_idx < trajectories.length; system_idx++) {
            System.out.println(String.format("System: %s", TimeSeriesGenerator.system_names[system_idx]));
            double[][][] trajectory1 = trajectories[system_idx];
            for (double[][] trajectory : trajectory1) {
                DRQA drqa = new DRQA(trajectory, trajectory, 0.05);
                drqa.computeRQA(2, 2, 2);
//                rqa_measures[system_idx * numberOfTimeSeriesPerClass + i] = new DoublePoint(drqa.allRQAMeasures(DRQA.HistogramStatistic.values()));
                rqa_measures.add(new DoublePoint(drqa.standardRQAMeasures()));
            }
        }
        System.out.println(String.format("Standard RQA of first point:\n%s", Arrays.toString(rqa_measures.get(0).getPoint())));

        // ------------------------------------
        // Cluster
        // ------------------------------------

        long before = System.currentTimeMillis();

        RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        final EuclideanDistance euclideanDistance = new EuclideanDistance();
        ClusterEvaluator<DoublePoint> evaluator = new SumOfClusterVariances<DoublePoint>(euclideanDistance);

        ArrayList<DoublePoint> sample = new ArrayList<DoublePoint>(sample_size);
        System.out.println(String.format("sample_size: %s", sample_size));

        // generate random subsamples of the data set
        for (int resample_idx = 0; resample_idx < numberOfResamples; resample_idx++) {

            System.out.println(String.format("Resample: %s", resample_idx+1));

            // generate random sample of half size
            Object[] randomSample = randomDataGenerator.nextSample(rqa_measures, sample_size);
            sample.clear(); for(Object point : randomSample) sample.add((DoublePoint) point);

            int bestK = 0;
            double bestScore = Double.POSITIVE_INFINITY;
            // cluster using different numbers of clusters
            for (int k = 2; k < 17; k++) {

                // seed the algorithm several times to have a better chance to escape local maxima
                // this can also be done using org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer,
                // but it doesn't return the actual scores, which is important in determining the number of clusters.
                double bestInitScore = Double.POSITIVE_INFINITY;
                for (int init_num = 0; init_num < numberOfInitializations; init_num++) {
                    KMeansPlusPlusClusterer<DoublePoint> clusterer = new KMeansPlusPlusClusterer<DoublePoint>(k, maxIterations, euclideanDistance, new ISAACRandom(init_num));
                    List<CentroidCluster<DoublePoint>> result = clusterer.cluster(sample);

                    final double score = evaluator.score(result);
                    bestInitScore = Math.min(bestInitScore, score);
                } // for initializations
                if(bestInitScore < bestScore){
                    bestScore = bestInitScore;
                    bestK = k;
                }
//                for(CentroidCluster<DoublePoint> cc : result) System.out.println(Arrays.toString(cc.getCenter().getPoint()));
            } // for number of clusters
            System.out.println(String.format("Best k=%s, score=%s", bestK, bestScore));
        } // for resampling

        System.out.println(String.format("Time for clustering: %s", System.currentTimeMillis() - before));


    }

    public void testLocate(){

    }

}