/*
 * Copyright (c) 2012 Oracle and/or its affiliates.
 * All rights reserved. Use is subject to license terms.
 *
 * This file is available and licensed under the following license:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the distribution.
 *  - Neither the name of Oracle nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package de.uberspace.wittcarl.experimental;

import de.uberspace.wittcarl.executable.BaseExecutable;
import javafx.animation.Animation;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Paint;
import javafx.stage.Stage;
import javafx.util.Duration;

public class Trajectory2DVisualization extends Application {

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    public void start(Stage primaryStage) {
        primaryStage.setTitle("Trajectory 2D");
        Group root = new Group();
        Canvas canvas = new Canvas(1020, 630);
        final GraphicsContext gc = canvas.getGraphicsContext2D();
        root.getChildren().add(canvas);
        primaryStage.setScene(new Scene(root));
        primaryStage.show();

        final double[][] trajectoryF = BaseExecutable.readFile("data/bern-barcelona/Data_F_Ind0002.txt", ",");
        final double[][] trajectoryN = BaseExecutable.readFile("data/bern-barcelona/Data_N_Ind0002.txt", ",");

        final Timeline tl = new Timeline();
        tl.setCycleCount(Animation.INDEFINITE);
        KeyFrame redraw = new KeyFrame(Duration.seconds(1. / 100),
                new EventHandler<ActionEvent>() {
                    int step = 1001;

                    public void handle(ActionEvent event) {
                        drawTrajectories(gc, trajectoryF, trajectoryN, step++);
                        if (step == trajectoryF[0].length) tl.stop();
                    }
                });
        tl.getKeyFrames().add(redraw);
        tl.play();


    }

    /**
     * Draws a series of basic shapes on the specified GraphicsContext.
     *
     * @param gc The GraphicsContext object to draw on
     */
    private void drawTrajectories(GraphicsContext gc, double[][] trajectoryF, double[][] trajectoryN, int step) {
        gc.clearRect(0, 0, 1e10, 1e10);
//        gc.setFill(Color.GREEN);
//        gc.setStroke(Color.BLUE);
//        gc.setLineWidth(5);
        if (step == 0) return;
        double offsetY = 300;
        double offsetX_F = 300;
        double offsetX_N = 600;
        gc.setStroke(Paint.valueOf("#ff0000ff"));
        gc.strokeOval(offsetX_F + trajectoryF[0][step], offsetY + trajectoryF[1][step], 5, 5);
        gc.strokeOval(offsetX_N + trajectoryN[0][step], offsetY + trajectoryN[1][step], 5, 5);

        gc.setStroke(Paint.valueOf("#00000033"));
        for (int i = 0; i < 1000; i++) {
            gc.strokeLine(offsetX_F + trajectoryF[0][step - i], offsetY + trajectoryF[1][step - i], offsetX_F + trajectoryF[0][step - i - 1], offsetY + trajectoryF[1][step - i - 1]);
            gc.strokeLine(offsetX_N + trajectoryN[0][step - i], offsetY + trajectoryN[1][step - i], offsetX_N + trajectoryN[0][step - i - 1], offsetY + trajectoryN[1][step - i - 1]);
        }

    }
}
