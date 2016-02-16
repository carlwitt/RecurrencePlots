package de.uberspace.wittcarl.gui;

import de.uberspace.wittcarl.Cacheable;
import de.uberspace.wittcarl.DRQA;
import de.uberspace.wittcarl.datagenerator.TimeSeriesGenerator;
import de.uberspace.wittcarl.phasespace.PhaseSpaceDistribution;
import de.uberspace.wittcarl.phasespace.PhaseSpaceReconstructed;
import de.uberspace.wittcarl.executable.BaseExecutable;
import javafx.beans.binding.Bindings;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.embed.swing.SwingFXUtils;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.geometry.Rectangle2D;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.*;
import javafx.scene.image.ImageView;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyEvent;
import javafx.scene.input.ZoomEvent;
import javafx.scene.layout.Pane;
import javafx.stage.FileChooser;
import javafx.stage.Screen;
import javafx.stage.Stage;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.awt.image.BufferedImage;
import java.io.File;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.ResourceBundle;

public class RPWindowController implements Initializable {

    Stage stage;
    FileChooser fileChooser = new FileChooser();

    double[][] trajectory;
    /** A version of the trajectory that accounts for range selection and embedding. */
    Cacheable<double[][]> transformedTrajectory;
    DRQA drqa;

    private final DoubleProperty recurrenceThreshold = new SimpleDoubleProperty();
    private final IntegerProperty embeddingDimension = new SimpleIntegerProperty();
    private final IntegerProperty embeddingDelay = new SimpleIntegerProperty();
    private final IntegerProperty considerPointsUpTo = new SimpleIntegerProperty();
    private final IntegerProperty eachKthPoint = new SimpleIntegerProperty(1);

    private final DoubleProperty noiseRatio = new SimpleDoubleProperty();

    /** Make it zoomable https://gist.github.com/james-d/7252698 */
    @FXML private LineChart<Double, Double> timeSeriesChart;
    @FXML private LineChart<Double, Double> subsequentDistanceDistributionChart;
    @FXML private LineChart<Double, Double> lineLengthHistogram;

    @FXML private ChoiceBox lineLengthTypeSelector;
    @FXML private ChoiceBox distanceDistributionSelector;

    @FXML private Pane mainWindowRoot;
    @FXML private ProgressBar progressBar;
    @FXML private Label progressLabel;
    @FXML private Label inputLengthLabel;
    @FXML private Slider recurrenceThresholdSlider;
    @FXML private Slider subsampleLengthSlider;
    @FXML private TextField recurrenceThresholdTextField;
    @FXML private TextField embeddingDimensionTextField;
    @FXML private TextField embeddingDelayTextField;
    @FXML private TextField crtLimit;
    @FXML private CheckBox sliderAutoUpdate;
    @FXML private CheckBox useDRQA;
    @FXML private Button computeRPButton;
    @FXML private CheckBox logScaleCheckBox;
    @FXML private TabPane crtTabPane;
    @FXML private TextArea crtStats;
    @FXML private TextArea rqaMeasures;
    @FXML private Slider noiseSlider;
    @FXML private TextField eachKthPointTextField;

    /** The zoom factor that is applied to the image. */
    private double imageScale = 1.0;
    @FXML private ImageView rpImageView;
    @FXML private ImageView crtImageView;


    /**
     * Called after the controls have been parsed from the XML. Sets up logic and components that could not be set up using the GUI builder.
     */
    @Override public void initialize(URL url, ResourceBundle rb) {

        // the transformedTrajectory provides a view of the trajectory that reflects the selected transformedTrajectory
        transformedTrajectory = new Cacheable<double[][]>() {
            private int eachKthPoint;
            private int embeddingDimension, embeddingDelay;
            private double noiseRatio;

            private int getSubsampleLength(){
                int prefixLength = (int) (subsampleLengthSlider.getValue() * trajectory[0].length);
                // per started block of k elements, one output element will be generated
                int k = getEachKthPoint();
                // full blocks + 1 block if there is a fractional block
                return prefixLength/k + ( prefixLength % k > 0 ? 1 : 0 );
            }
            @Override public synchronized boolean isValid() {
                if(cachedValue == null) return false;
                if( eachKthPoint != getEachKthPoint()) return false;
                // any embedding dimension <= 0 signifies that no embedding should be used, thus the actual dimension does not matter
                if( getEmbeddingDimension() > 0 && (embeddingDimension != getEmbeddingDimension())) return false;
                // for an embedding dimension of <= 1, the delay is insignificant
                if( getEmbeddingDimension() > 1 && (embeddingDelay != getEmbeddingDelay())) return false;
                if( getNoiseRatio() != noiseRatio) return false;
                return trajectory == null || cachedValue == null || getSubsampleLength() == cachedValue[0].length;
            }
            @Override public synchronized void recompute() {
                if(trajectory == null) cachedValue = null;
                else{
                    // create an empty array with the desired number of dimensions
                    cachedValue = new double[trajectory.length][];
                    eachKthPoint = getEachKthPoint();

                    // crop to sampling size
                    int newLength = getSubsampleLength();
                    for (int dimIdx = 0; dimIdx < trajectory.length; dimIdx++) {
                        cachedValue[dimIdx] = new double[newLength]; // Arrays.copyOf(trajectory[dimIdx], newLength);
                        for (int i=0, t=0; i < newLength; i++, t+=eachKthPoint) {
                            cachedValue[dimIdx][i] = trajectory[dimIdx][t];
                        }
                    }

                    int dim = getEmbeddingDimension();
                    int tau = getEmbeddingDelay();
                    double noiseRatio = getNoiseRatio();

                    if(dim > 0 && tau > 0){
                        cachedValue = PhaseSpaceReconstructed.embed(cachedValue[0], dim, tau);
                    }

                    cachedValue = TimeSeriesGenerator.addNoise(cachedValue, 0.05, noiseRatio);

                    this.noiseRatio = noiseRatio;
                    this.embeddingDimension = dim;
                    this.embeddingDelay = tau;
                }
            }
        };

        // sync recurrence threshold slider and text field
        Bindings.bindBidirectional(recurrenceThresholdTextField.textProperty(), recurrenceThresholdSlider.valueProperty(), DecimalFormat.getInstance());
        // sync GUI and model recurrence threshold
        Bindings.bindBidirectional(recurrenceThresholdTextField.textProperty(), recurrenceThresholdProperty(), DecimalFormat.getInstance());

        // sync GUI and model embedding parameters
        Bindings.bindBidirectional(embeddingDimensionTextField.textProperty(), embeddingDimensionProperty(), DecimalFormat.getIntegerInstance());
        Bindings.bindBidirectional(embeddingDelayTextField.textProperty(), embeddingDelayProperty(), DecimalFormat.getIntegerInstance());

        // sync GUI and noise parameter
        Bindings.bindBidirectional(noiseSlider.valueProperty(), noiseRatioProperty());

        Bindings.bindBidirectional(eachKthPointTextField.textProperty(), eachKthPointProperty(), DecimalFormat.getIntegerInstance());

        // enable the compute button only if the auto-update checkbox is not on.
        computeRPButton.disableProperty().bind(sliderAutoUpdate.selectedProperty());

        // recompute RP on parameter changes
        subsampleLengthSlider.valueProperty().addListener(this::parametersChanged);
        eachKthPointTextField.textProperty().addListener((obs, ov, nv) -> { parametersChanged(null, null, null); });
        recurrenceThresholdProperty().addListener(this::parametersChanged);
        embeddingDimensionProperty().addListener(this::parametersChanged);
        embeddingDelayProperty().addListener(this::parametersChanged);
        noiseRatioProperty().addListener(this::parametersChanged);

        // Make CRT controls update the computation
        // size of the CRT histogram
        crtLimit.textProperty().addListener(new ChangeListener<String>() {
            @Override public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
                parametersChanged(null, Integer.parseInt(oldValue), Integer.parseInt(newValue));
            }
        });
        // CRT log scale option
        logScaleCheckBox.selectedProperty().addListener(new ChangeListener<Boolean>() {
            @Override public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                if(oldValue != newValue) parametersChanged(null, 0 , 0);
            }
        });

        // make the CRT image use all the available vertical space
        crtImageView.fitWidthProperty().bind(crtTabPane.widthProperty().subtract(20));
        crtImageView.fitHeightProperty().bind(crtTabPane.heightProperty());

        // swap the data of the line length histogram upon selecting another line type
        lineLengthTypeSelector.setItems(FXCollections.observableArrayList("DIAGONAL", "VERTICAL", "WHITE_VERTICAL", "ORTHOGONAL"));//DRQA.LineType.DIAGONAL, DRQA.LineType.VERTICAL, DRQA.LineType.WHITE_VERTICAL, DRQA.LineType.ORTHOGONAL
        lineLengthTypeSelector.getSelectionModel().selectedIndexProperty().addListener(this::updateLineLengthHistogram);
        lineLengthTypeSelector.getSelectionModel().select(0);

        distanceDistributionSelector.setItems(FXCollections.observableArrayList("SUBSEQ", "PAIRWISE"));
        distanceDistributionSelector.getSelectionModel().selectedIndexProperty().addListener((obs, ov, nv)-> updateDistanceDistributionChart());
        distanceDistributionSelector.getSelectionModel().select(0);

        useDRQA.selectedProperty().addListener((obs,ov,nv)->updateLineLengthHistogram(null, null, lineLengthTypeSelector.getSelectionModel().getSelectedIndex()));


    }

    public double[][] getTransformedTrajectory(){
        return transformedTrajectory.get();
    }

    /**
     * Called whenever the recurrence threshold or embedding parameters change.
     */
    @FXML void parametersChanged( ObservableValue<? extends Number> observable, Number oldValue, Number newValue ){
        if(sliderAutoUpdate.isSelected()) this.compute();
    }

    @FXML public void compute(){
        double[][] transformedTrajectory = getTransformedTrajectory();
        if(transformedTrajectory == null) return;

        // compute and display RP
        BufferedImage rp = DRQA.getRPImage(transformedTrajectory, transformedTrajectory, recurrenceThresholdSlider.getValue());
        rpImageView.setImage(SwingFXUtils.toFXImage(rp, null));
        applyImageScale();

        // compute and display CRT
        DRQA.conditional_ww_limit = Integer.parseInt(crtLimit.getText());
        DRQA.CRT_LOG_SCALE = logScaleCheckBox.isSelected();
        drqa = new DRQA(transformedTrajectory, transformedTrajectory, recurrenceThresholdSlider.getValue());
        BufferedImage crt = drqa.getCRTImage(DRQA.conditional_ww_limit, drqa.conditional_ww);
        crtImageView.setImage(SwingFXUtils.toFXImage(crt, null));
        String[] stats = drqa.crtStatistics().split("\t");
        crtStats.setText(String.format("mean row: %.2f\tmean col: %.2f\ncorrelation: %.2f\nmax row: %s\tmax col: %s\nlocal maxima: %s\nentropy: %.2f",Double.parseDouble(stats[0]), Double.parseDouble(stats[1]), Double.parseDouble(stats[2]),stats[3],stats[4],stats[5],Double.parseDouble(stats[6])));

        drqa.computeRQA(2,2,2);
        rqaMeasures.setText(drqa.printableString(DRQA.STANDARD_RQA));

        updateTimeSeriesChart();
        updateDistanceDistributionChart();
        updateLineLengthHistogram(null, null, lineLengthTypeSelector.getSelectionModel().getSelectedIndex());
    }

    @FXML
    public void openFile() {
        stage = (Stage) mainWindowRoot.getScene().getWindow();
        fileChooser.setSelectedExtensionFilter(new FileChooser.ExtensionFilter("Comma separated text files", "*.txt, *.csv"));
        File in = fileChooser.showOpenDialog(stage);
        if (in != null) {
            String filename = in.getPath();
            setTrajectory(BaseExecutable.readFile(filename, ","));
        }
    }

    public void setTrajectory(double[][] trajectory) {
        this.trajectory = trajectory;
        inputLengthLabel.setText("Input Length: " + trajectory[0].length);
        compute();
    }

    @FXML public void saveRPImage(){
        throw new NotImplementedException();
//        fileChooser.setSelectedExtensionFilter(new FileChooser.ExtensionFilter("PNG Image (*.png)", "*.png"));
//        File out = fileChooser.showSaveDialog(stage);
        // if(out != null) out.getPath();
    }

    @FXML public void saveCRTImage(){
        throw new NotImplementedException();
//        fileChooser.setSelectedExtensionFilter(new FileChooser.ExtensionFilter("PNG Image (*.png)", "*.png"));
//        File out = fileChooser.showSaveDialog(stage);
        // if(out != null) out.getPath();
    }

    public void showWindow() {
        // resize primary stage to full screen
        Screen primaryScreen = Screen.getPrimary();
        Rectangle2D bounds = primaryScreen.getVisualBounds();
        Stage mainWindowStage = (Stage) mainWindowRoot.getScene().getWindow();
        mainWindowStage.setX(bounds.getMinX());
        mainWindowStage.setY(bounds.getMinY()+98);
        mainWindowStage.setWidth(bounds.getWidth());
        mainWindowStage.setHeight(1080);//bounds.getHeight()
        mainWindowStage.setOnCloseRequest(event -> quit());
    }

    private void updateDistanceDistributionChart() {

        double[][] transformedTrajectory = getTransformedTrajectory();
        if(transformedTrajectory == null) return;

        double approxDiameter = PhaseSpaceDistribution.maxPhaseSpaceDiameterApproximate(transformedTrajectory);

        long[] histo;
        // subsequent distance distribution
        if(distanceDistributionSelector.getSelectionModel().getSelectedIndex() == 0){
            histo = PhaseSpaceDistribution.approximateSubsequentDistanceDistribution(transformedTrajectory, 1500, approxDiameter);
        } else {
            histo = PhaseSpaceDistribution.approximateDistanceDistribution(transformedTrajectory, 1500, approxDiameter);
        }


        subsequentDistanceDistributionChart.getData().clear();
        XYChart.Series series = new XYChart.Series();
        for (int i = 0; i < histo.length; i++) {
            series.getData().add(new XYChart.Data(i, histo[i]));
        }
        subsequentDistanceDistributionChart.getData().add(series);

    }

    private void updateTimeSeriesChart() {
        // clear time series chart
        timeSeriesChart.getData().clear();

        // fill time series chart
        double[][] transformedTrajectory = getTransformedTrajectory();
        for (int i = 0; i < transformedTrajectory.length; i++) {
            XYChart.Series series = new XYChart.Series();
            series.setName("Dimension "+(i+1));
            for (int j = 0; j < transformedTrajectory[i].length; j++) {
                series.getData().add(new XYChart.Data(j, transformedTrajectory[i][j]));
            }
            timeSeriesChart.getData().add(series);
        }
    }

    private void updateLineLengthHistogram(ObservableValue obs, Number oldValue, Number newValue){
        if(drqa == null) return;

        long[][] histos;
        if(useDRQA.isSelected())
            histos = new long[][]{drqa.l_hist, drqa.v_hist, drqa.w_hist, drqa.r_hist};
        else{
            // add indefinite and definite histograms
            histos = new long[][]{Arrays.copyOf(drqa.l_hist,drqa.l_hist.length), Arrays.copyOf(drqa.v_hist,drqa.v_hist.length), Arrays.copyOf(drqa.w_hist,drqa.w_hist.length), Arrays.copyOf(drqa.r_hist,drqa.r_hist.length),};
            long[][] indefiniteHistos = new long[][]{drqa.l_hist_indefinite, drqa.v_hist_indefinite, drqa.w_hist_indefinite, drqa.r_hist_indefinite};
            // add the indefinite histos to the definite histos to gain traditional rqa
            for (int i = 0; i < histos.length; i++) {
                for (int j = 0; j < indefiniteHistos[i].length; j++) histos[i][j] += indefiniteHistos[i][j];
            }
        }

        lineLengthHistogram.getData().clear();
        long[] histogram = histos[newValue.intValue()];

        // fill time series chart
        XYChart.Series series = new XYChart.Series();
        for (int i = 0; i < histogram.length; i++) {
            series.getData().add(new XYChart.Data(i, histogram[i]));
        }
        lineLengthHistogram.getData().add(series);
    }

    private void applyImageScale(){
        rpImageView.setFitWidth(rpImageView.getImage().getWidth()*imageScale);
        rpImageView.setFitHeight(rpImageView.getImage().getHeight()*imageScale);
    }
    @FXML public void resetImageZoom(){
        imageScale = 1;
        applyImageScale();
    }
    @FXML public void imageZoom(ZoomEvent event){
        event.consume();
        imageScale *= event.getZoomFactor();
        applyImageScale();
    }

    /** Allow for arrow key increase/decrease of the value of text fields that store integers. */
    @FXML void keyArrowsTextManipulator(KeyEvent event){
        if(event.getSource() instanceof TextField){
            TextField tf = (TextField) event.getSource();
            int oldValue = Integer.parseInt(tf.getText());
            if(event.getCode() == KeyCode.DOWN) tf.setText("" + (oldValue-1));
            if(event.getCode() == KeyCode.UP) tf.setText("" + (oldValue+1));
        }
    }

    public int getEmbeddingDimension() {
        return embeddingDimension.get();
    }

    public IntegerProperty embeddingDimensionProperty() {
        return embeddingDimension;
    }

    public void setEmbeddingDimension(int embeddingDimension) {
        this.embeddingDimension.set(embeddingDimension);
    }
    public int getEmbeddingDelay() {
        return embeddingDelay.get();
    }

    public IntegerProperty embeddingDelayProperty() {
        return embeddingDelay;
    }

    public void setEmbeddingDelay(int embeddingDelay) {
        this.embeddingDelay.set(embeddingDelay);
    }

    public double getRecurrenceThreshold() {
        return recurrenceThreshold.get();
    }

    public DoubleProperty recurrenceThresholdProperty() {
        return recurrenceThreshold;
    }

    public void setRecurrenceThreshold(double recurrenceThreshold) {
        this.recurrenceThreshold.set(recurrenceThreshold);
    }

    public int getConsiderPointsUpTo() {
        return considerPointsUpTo.get();
    }

    public IntegerProperty considerPointsUpToProperty() {
        return considerPointsUpTo;
    }

    public void setConsiderPointsUpTo(int considerPointsUpTo) {
        this.considerPointsUpTo.set(considerPointsUpTo);
    }


    public double getNoiseRatio() { return noiseRatio.get(); }

    public DoubleProperty noiseRatioProperty() { return noiseRatio; }

    public void setNoiseRatio(double noiseRatio) { this.noiseRatio.set(noiseRatio); }
    public int getEachKthPoint() { return eachKthPoint.get(); }
    public IntegerProperty eachKthPointProperty() { return eachKthPoint; }
    public void setEachKthPoint(int eachKthPoint) { this.eachKthPoint.set(eachKthPoint); }

    public void quit(){
        System.exit(0);
    }

}
